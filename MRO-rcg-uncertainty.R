#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there are exactly three arguments: if not, return an error
if (length(args)!=4) {
  stop("Need four arguments: a Phyllip or Newick tree, binary strings of loss patterns, an output folder, and an output label (1: 1b --> 0; 2: 1b --> 1).\n Usage: MRO-rcg-uncertainty.R [tree file] [barcodes] [output folder] [output label]  \n", call.=FALSE)
}

library(phytools)
library(parallel)
library(ggplot2)

treefile <- as.character(args[1])
datafile <- as.character(args[2])
output <- as.character(args[3])
tag <- as.character(args[4])

tree <- read.tree(paste("Data/", treefile,sep =""))
df <- read.csv(paste("Data/RandomizedData-2/", datafile, sep = ""))
outfolder <- paste("Results/", output,"/", sep = "")

# IGJ March 2025 -- updated this file
source("hypertraps.R")

# harmonise labels across tree and barcode data
# IGJ March 2025 -- think I added something here
# df$Organism = gsub("_", " ", df$Organism)

#tree$tip.label = gsub("'", "", tree$tip.label)
tree$tip.label = gsub("_", " ", tree$tip.label)

# wrapper function for HyperTraPS analysis
parallel.fn = function(seed, src.data, featurenames) {
  return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                    losses = 1, length = 5, kernel = 6,
                    seed = seed, samplegap = 100, penalty = 0.,
                    featurenames = featurenames))
}

# Set a seed
set.seed(1234)

# ancestral state reconstruction and transition gathering
data.ncbi = curate.tree(tree, df, losses=TRUE)

# run these experiments in parallel. should take a few core minutes each
n.seed = 3
featurenames = colnames(df)[2:length(colnames(df))]
parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                              MoreArgs = list(src.data=data.ncbi,
                                              featurenames=featurenames),
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), n.seed))

g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
g.blurb = plotHypercube.summary(parallelised.runs[[1]]) + theme(legend.position="none")
g.bub.1 = plotHypercube.bubbles(parallelised.runs[[1]]) + theme(legend.position="none")
g.bub.2 = plotHypercube.bubbles(parallelised.runs[[2]]) + theme(legend.position="none")
g.bub.3 = plotHypercube.bubbles(parallelised.runs[[3]]) + theme(legend.position="none")
g.sg.1 = plotHypercube.sampledgraph2(parallelised.runs[[1]], 
                                     thresh=0.05, use.arc=FALSE, 
                                     node.labels = FALSE,
                                     edge.label.size=4, no.times=TRUE) + 
  theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))
 
# compare outputs and check convergence
sf = 2
png(paste(outfolder, "mro-convergence-", tag, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
dev.off()

# check summary
sf = 2
png(paste(outfolder, "mro-trace-", tag, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
print(g.blurb)
dev.off()

# dynamics
#g.sg.1
# summary plot
g.summary = ggarrange(g.data, ggarrange(g.bub.1, g.sg.1, ncol=1, 
                                        heights=c(1,2), labels=c("B", "C")), 
                      nrow=1, widths=c(1.25,1), labels=c("A",""))
#g.summary

sf = 2
png(paste(outfolder, "mro-summary-", tag, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
print(g.summary)
dev.off()

# Check influences
g.influences <- plotHypercube.influences(parallelised.runs[[1]])
g.influencegraph <- plotHypercube.influencegraph(parallelised.runs[[1]])

g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))

sf = 2
png(paste(outfolder, "mro-influences-", tag, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
print(g.infl)
dev.off()

r = parallelised.runs[[1]]$routes
r.ci = r[r[,1]==0,]
r.tca = r[r[,1]==7,]
r.c34 = r[r[,1]==2 | r[,1]==3,]
r.pdh = r[r[,1]==5,]
r.other = r[r[,1]!=0,]

#save.image(paste(outfolder, "mro-save-", tag, ".RData",sep = ""))

tmp.ci = tmp.other = tmp.tca = tmp.c34 = tmp.pdh = parallelised.runs[[1]]
tmp.ci$routes = r.ci
tmp.tca$routes = r.tca
tmp.c34$routes = r.c34
tmp.pdh$routes = r.pdh
tmp.other$routes = r.other
g.pathways = ggarrange(plotHypercube.motifs(parallelised.runs[[1]])+theme(legend.position="none"),
          plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.c34)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.pdh)+theme(legend.position="none"),
          ncol=1, labels=c("A", "B", "C", "D", "E"))

sf = 2
png(paste(outfolder, "mro-pathways-", tag, ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
print(g.pathways)
dev.off()

