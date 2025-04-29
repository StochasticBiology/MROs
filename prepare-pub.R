#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there are exactly two arguments: if not, return an error
if (length(args)!=2) {
  stop("Need two arguments: the integer label of the results/input folder and the cv threshold for the summary plot.\n Usage: MRO-rcg.R [integer label (1: 1b --> 0; 2: 1b --> 1)] [cv threshold]\n", call.=FALSE)
}


library(igraph)
library(ggraph)
library(ggtree)
library(ggtreeExtra)
library(ggpubr)
source("hypertraps.R")

mytag = as.character(args[1])
cv.thresh = as.numeric(args[2])

#### Figure 1

load(paste("Results/HyperTraPS-", mytag, "/mro-save-0-0.RData",sep=""))
inter.set = parallelised.runs[[1]]
igj.plot.2 = plotHypercube.influencegraph(inter.set, cv.thresh = cv.thresh, label.size = 4)

igj.plot.2a = plotHypercube.curated.tree(data.ncbi, names=TRUE) 

igj.plot.2b = plotHypercube.bubbles(inter.set) + theme(legend.position="none")
igj.plot.2c = plotHypercube.sampledgraph2(inter.set, node.labels = FALSE,
                                          no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")

ggarrange(igj.plot.2a +
            scale_x_continuous(expand = c(0.5,0)), 
          ggarrange(igj.plot.2b, igj.plot.2c, nrow=2, heights=c(1,1.4)), 
          nrow=1, widths=c(2,1))

summary.plot = ggarrange(igj.plot.2a +
            scale_x_continuous(expand = expansion(mult = c(0.6, 0), 
                   add = c(-2, 0))), 
                   ggarrange(igj.plot.2b, igj.plot.2c, 
                             nrow=2, heights=c(1,1.4),
                             labels=c("B", "C")), 
                   nrow=1, widths=c(1.,1), labels=c("A", ""))

sf = 2
png("Results/fig-1.png", width=800*sf, height=600*sf, res=72*sf)
print(summary.plot)
dev.off()

### Fig 2: Influences between features, and evolutionary pathway structure
r = parallelised.runs[[1]]$routes
r.c34tca = r[r[,1] == 2 | r[,1] == 3 | r[,1] == 7,]
r.c1 = r[r[,1]==0,]

tmp.c1 = tmp.c34tca = parallelised.runs[[1]]
tmp.c34tca$routes = r.c34tca
tmp.c1$routes = r.c1
g.pathways = ggarrange(plotHypercube.motifs(parallelised.runs[[1]])+theme(legend.position="none"),
                       plotHypercube.motifs(tmp.c1)+theme(legend.position="none"),
                       plotHypercube.motifs(tmp.c34tca)+theme(legend.position="none"),
                       ncol=1, labels=c("i", "ii", "iii"))

g.influences = plotHypercube.influencegraph(inter.set, cv.thresh = cv.thresh, label.size = 4)

#### 3 AM: can't get these labels right; might want to separate
fig.2 = ggarrange(g.pathways, g.influences, labels = c("A","B"), ncol = 2)

sf = 2
png("Results/fig-2.png", width=800*sf, height=600*sf, res=72*sf)
print(fig.2)
dev.off()


#### Fig. 3: ciliates vs. apicomplexans
load(paste("Results/Ciliophora-", mytag, "/mro-save-0-0.RData",sep=""))
set.1 = parallelised.runs[[1]]
load(paste("Results/Apicomplexans-", mytag ,"/mro-save-0-0.RData",sep=""))
set.2 = parallelised.runs[[1]]

fig.3a.1 = plotHypercube.bubbles(set.1) + theme(legend.position="none")
fig.3a.2 = plotHypercube.sampledgraph2(set.1, node.labels = FALSE,
                                          no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")
fig.3b.1 = plotHypercube.bubbles(set.2) + theme(legend.position="none")
fig.3b.2 = plotHypercube.sampledgraph2(set.2, node.labels = FALSE,
                                          no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")

#### 3 AM: can't get these labels right
fig.3 = ggarrange(ggarrange(fig.3a.1,fig.3a.2, nrow = 2, labels = c("i","ii")) + theme(legend.position = "none"),
                  ggarrange(fig.3b.1,fig.3b.2, nrow = 2, labels = c("i","ii")) + theme(legend.position = "none"),
                  nrow = 1, labels = c("A","B"))
sf = 2
png("Results/fig-3.png", width=800*sf, height=600*sf, res=72*sf)
print(fig.3)
dev.off()

#### SI Fig 1: HyperTraPS-CT inference
load(paste("Results/HyperTraPS-CT-", mytag, "/mro-save-0-1.RData",sep = ""))

inter.set = parallelised.runs[[1]]
igj.plot.2 = plotHypercube.influencegraph(inter.set, cv.thresh = cv.thresh, label.size = 4)

igj.plot.2a = plotHypercube.curated.tree(data.ncbi, names=TRUE) 

igj.plot.2b = plotHypercube.bubbles(inter.set) + theme(legend.position="none")
igj.plot.2c = plotHypercube.sampledgraph2(inter.set, node.labels = FALSE,
                                          no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")

ggarrange(igj.plot.2a +
            scale_x_continuous(expand = c(0.5,0)), 
          ggarrange(igj.plot.2b, igj.plot.2c, nrow=2, heights=c(1,1.4)), 
          nrow=1, widths=c(2,1))

summary.plot = ggarrange(igj.plot.2a +
            scale_x_continuous(expand = expansion(mult = c(0.6, 0), 
                   add = c(-2, 0))), 
                   ggarrange(igj.plot.2b, igj.plot.2c, 
                             nrow=2, heights=c(1,1.4),
                             labels=c("B", "C")), 
                   nrow=1, widths=c(1.,1), labels=c("A", ""))

sf = 2
png(paste("Results/SI-fig-1.png", sep = ""), width=800*sf, height=600*sf, res=72*sf)
print(summary.plot)
dev.off()

#### SI fig. 2: Robustness wrt uncertainty in source data
# pull all 10 resampled inference runs into a single list
unc.list = list()
for(iains.i in 1:10) {
  load(paste("Results/HyperTraPS-uncertainty-", mytag, "/mro-save-", iains.i, "-0.RData", sep=""))
  unc.list[[iains.i]] = parallelised.runs[[1]]
}

# construct a comparison plot for this list
raw.plot = plotHypercube.bubbles.compare(unc.list) 
f.names = unc.list[[1]]$featurenames
png("Results/SI-fig-2.png", width = 800*sf, height = 600*sf, res = 72*sf)
raw.plot + scale_fill_viridis_d(option="magma") + scale_y_continuous(breaks = 0:(length(f.names)-1), labels = f.names)
dev.off()

# SI fig. 3: Iain


# SI fig. 4: 
load(paste("Results/HyperTraPS-CT-", mytag,"/mro-save-0-1.RData", sep = ""))
set.1 = which(parallelised.runs[[1]]$routes[,1] == 0)
set.2 = which(parallelised.runs[[1]]$routes[,1] %in% c(2,3,7))
times.df = rbind(
     data.frame(pathset = "CI", 
                step = 1:ncol(parallelised.runs[[1]]$times),
                means = apply(log10(parallelised.runs[[1]]$times[set.1,]), 2, mean), 
                sds = apply(log10(parallelised.runs[[1]]$times[set.1,]), 2, sd)),
     data.frame(pathset = "alt", 
                step = 1:ncol(parallelised.runs[[1]]$times),
                means = apply(log10(parallelised.runs[[1]]$times[set.2,]), 2, mean), 
                sds = apply(log10(parallelised.runs[[1]]$times[set.2,]), 2, sd)) 
 )
timingplot <- ggplot(times.df, aes(x=step, ymin=means-sds, ymax=means+sds, color=pathset)) + 
  geom_errorbar(width=0.2, position="dodge")+
  labs(x = "Step", y = "log10 time (mean +- sd)", color = "Pathway") +
  theme_minimal()+scale_x_continuous(breaks=seq(-1,10,by=1))

sf = 2
png("Results/SI-fig-4.png", width=300*sf,height=200*sf,res=72*sf)
timingplot
dev.off()
