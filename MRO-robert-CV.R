#!/usr/bin/env Rscript

#### Script to crossvalidate the findings from original HyperTraPS (note: not CT)
<<<<<<< Updated upstream
=======
args = commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Need two arguments, a Phyllip or Newick tree and barcodes of loss patterns.\n", call.=FALSE)
}

>>>>>>> Stashed changes
library(phytools)
library(parallel)
library(sybilSBML)
library(ggplot2)

<<<<<<< Updated upstream
#### Run HyperTraPS without timings on the NCBI tree 10 times with 10% of loss patterns associated with observations randomly changed ####
##   Save each of the data sets and all of the runs in RData files labelled according to the index of validation sets

#### Run HyperTraPS without timings on the NCBI tree ####
source("hypertraps.R")

# read data
tree = read.tree("Data/mro-ncbi-tree-2025.phy")
df = read.csv("Data/mro-barcodes-2025.csv")

# harmonise labels across tree and barcode data
# IGJ March 2025 -- think I added something here
df$Organism = gsub("_", " ", df$Organism)

tree$tip.label = gsub("'", "", tree$tip.label)
tree$tip.label = gsub("_", " ", tree$tip.label)

# ancestral state reconstruction and transition gathering
data.ncbi = curate.tree(tree, df, losses=TRUE)

# wrapper function for HyperTraPS analysis
parallel.fn = function(seed, src.data, featurenames) {
  return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                    losses = 1, length = 5, kernel = 6,
                    seed = seed, samplegap = 10, penalty = 0.,
                    featurenames = featurenames))
}

# run these experiments in parallel. should take a few core minutes each
n.seed = 5
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
png("mro-CV-outputs.png", width = 700*sf, height = 600*sf, res = 72*sf)
print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
dev.off()

# check summary
sf = 2
png("mro-CV-convergence.png", width = 700*sf, height = 600*sf, res = 72*sf)
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
png("mro-summary.png", width=700*sf, height=600*sf, res=72*sf)
print(g.summary)
dev.off()

r = parallelised.runs[[1]]$routes
r.ci = r[r[,1]==0,]
r.tca = r[r[,1]==7,]
r.c34 = r[r[,1]==2 | r[,1]==3,]
r.other = r[r[,1]!=0,]

tmp.ci = tmp.other = tmp.tca = tmp.c34 = parallelised.runs[[1]]
tmp.ci$routes = r.ci
tmp.tca$routes = r.tca
tmp.c34$routes = r.c34
tmp.other$routes = r.other
g.pathways = ggarrange(plotHypercube.motifs(parallelised.runs[[1]])+theme(legend.position="none"),
          plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.c34)+theme(legend.position="none"), 
          ncol=1, labels=c("A", "B", "C", "D"))

sf = 2
png("mro-CV-pathways.png", width=600*sf, height=600*sf, res=72*sf)
print(g.pathways)
dev.off()
=======
source("hypertraps.R")

# Read data
treefile <- as.character(args[1])
datafile <- as.character(args[2])

tree <- read.tree(paste("Data/",treefile,sep =""))
df <- read.csv(paste("Data/",datafile, sep = ""))

#### Run HyperTraPS without timings on the NCBI tree 10 times with 10% of loss patterns associated with observations randomly changed ####
##   Save each of the data sets and all of the runs in RData files labelled according to the index of validation sets

# Set a seed
set.seed(1234)

for(k in 1:10){
  this_df <- df

  # Change 10% of the observations to the opposite (0s to 1s and 1s to 0s)
  change <- sample(1:nrow(this_df), size = round(nrow(this_df)/10))

  tmp_df <- this_df[change,]
  tmp_df[tmp_df == 0] = 1
  tmp_df[this_df[change,] == 1] = 0

  this_df[change,] <- tmp_df

  # Save this_df
  write.csv(this_df, file = paste("Results/HyperTraPS-uncertainty/randomized-dataset-",k,".csv",sep = ""))

  # harmonise labels across tree and barcode data
  # IGJ March 2025 -- think I added something here
  this_df$Organism = gsub("_", " ", this_df$Organism)

  tree$tip.label = gsub("'", "", tree$tip.label)
  tree$tip.label = gsub("_", " ", tree$tip.label)

  # ancestral state reconstruction and transition gathering
  data.ncbi = curate.tree(tree, this_df, losses=TRUE)

  # wrapper function for HyperTraPS analysis
  parallel.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                    losses = 1, length = 5, kernel = 6,
                    seed = seed, samplegap = 100, penalty = 0.,
                    featurenames = featurenames))
  }

  # run these experiments in parallel. should take a few core minutes each
  n.seed = 5
  featurenames = colnames(this_df)[2:length(colnames(this_df))]
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
  png(paste("Results/HyperTraPS-uncertainty/mro-CV-convergence-",k,".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
  print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
  dev.off()
  
  # check summary
  sf = 2
  png(paste("Results/HyperTraPS-uncertainty/mro-CV-trace-",k,".png",sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
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
  png(paste("Results/HyperTraPS-uncertainty/mro-CV-summary-",k,".png",sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.summary)
  dev.off()
  
  r = parallelised.runs[[1]]$routes
  r.ci = r[r[,1]==0,]
  r.tca = r[r[,1]==7,]
  r.c34 = r[r[,1]==2 | r[,1]==3,]
  r.other = r[r[,1]!=0,]
  
  tmp.ci = tmp.other = tmp.tca = tmp.c34 = parallelised.runs[[1]]
  tmp.ci$routes = r.ci
  tmp.tca$routes = r.tca
  tmp.c34$routes = r.c34
  tmp.other$routes = r.other
  g.pathways = ggarrange(plotHypercube.motifs(parallelised.runs[[1]])+theme(legend.position="none"),
            plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
            plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
            plotHypercube.motifs(tmp.c34)+theme(legend.position="none"), 
            ncol=1, labels=c("A", "B", "C", "D"))
  
  sf = 2
  png(paste("Results/HyperTraPS-uncertainty/mro-CV-pathways-",k,".png",sep= ""), width=600*sf, height=600*sf, res=72*sf)
  print(g.pathways)
  dev.off()

  k <- k + 1
}

save.image("Results/HyperTraPS-uncertainty/mro-save-CV.RData")
>>>>>>> Stashed changes
