library(phytools)
library(parallel)
source("hypertraps.R")

# read data
tree = read.tree("Data/mro-ncbi-tree.phy")
df = read.csv("Data/mro-barcodes.csv")

# harmonise labels across tree and barcode data
df$Organism = gsub("_", " ", df$Organism)
tree$tip.label = gsub("'", "", tree$tip.label)

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
n.seed = 3
featurenames = colnames(df)[2:length(colnames(df))]
parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                               MoreArgs = list(src.data=data.ncbi,
                                               featurenames=featurenames),
                               SIMPLIFY = FALSE,
                               mc.cores = min(detectCores(), n.seed))

g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
g.bub.1 = plotHypercube.bubbles(parallelised.runs[[1]]) + theme(legend.position="none")
g.bub.2 = plotHypercube.bubbles(parallelised.runs[[2]]) + theme(legend.position="none")
g.bub.3 = plotHypercube.bubbles(parallelised.runs[[3]]) + theme(legend.position="none")
g.sg.1 = plotHypercube.sampledgraph2(parallelised.runs[[1]], 
                                     thresh=0.05, use.arc=FALSE, 
                                     node.labels = FALSE,
                                     edge.label.size=4, no.times=TRUE) + 
  theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))

# compare outputs and check convergence
ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1)
# dynamics
g.sg.1
# summary plot
g.summary = ggarrange(g.data, ggarrange(g.bub.1, g.sg.1, ncol=1, 
                                        heights=c(1,2), labels=c("B", "C")), 
                      nrow=1, widths=c(1.25,1), labels=c("A",""))
g.summary

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
png("mro-pathways.png", width=600*sf, height=600*sf, res=72*sf)
print(g.pathways)
dev.off()
