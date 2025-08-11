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
library(stringr)
library(hypertrapsct)
#source("hypertraps.R")

mytag = as.character(args[1])
cv.thresh = as.numeric(args[2])

plotHypercube.curated.tree.mro = function(tree.set, 
                                      scale.fn = geom_treescale(y=20, linesize=3, width =0.01),
                                      names = FALSE,
                                      line.color = "#000000",
                                      text.color = "#000000",
                                      font.size=4,
                                      hjust=0) {
  data.m = tree.set$data[,2:ncol(tree.set$data)]
  rownames(data.m) = tree.set$data[,1]
  data.m = tree.set$data[1:length(tree.set$tree$tip.label), 2:ncol(tree.set$data)]
  rownames(data.m) = tree.set$data$label[1:length(tree.set$tree$tip.label)]
  data.m[which(data.m=="?", arr.ind = TRUE)] = 0.5
  data.m = apply(data.m, c(1,2), as.numeric)
  
  g.core = ggtree(tree.set$tree, color=line.color) + scale.fn
  if(names == TRUE) {
    g.core = ggtree(tree.set$tree, color=line.color, branch.length = "none") + scale.fn + 
      geom_tiplab(size=2.5, hjust=0, color=text.color)
  } else {
    ggtree(tree.set$tree, color=line.color) + scale.fn
  }
  this.plot = gheatmap(g.core, data.m, low="white", high="#AAAAAA",
                       colnames_angle=90, hjust=hjust, font.size=font.size, offset=5) +
    theme(legend.position="none")
  
  return(this.plot)
}

#### Figure 1

load(paste("Results/HyperTraPS-", mytag, "/mro-save-0-0.RData",sep=""))
inter.set = parallelised.runs[[1]]
igj.plot.2 = plotHypercube.influencegraph(inter.set, cv.thresh = cv.thresh, label.size = 4)

igj.plot.2a = plotHypercube.curated.tree.mro(data.ncbi, names=TRUE, 
                                         line.color="#88888888", text.color = "#000000",
                                         scale.fn = NULL) 

igj.plot.2b = plotHypercube.bubbles(inter.set) + theme(legend.position="none")
igj.plot.2c = plotHypercube.sampledgraph2(inter.set, node.labels = FALSE,
                                          no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")

ggarrange(igj.plot.2a +
            scale_x_continuous(expand = c(0.5,0)), 
          ggarrange(igj.plot.2b, igj.plot.2c, nrow=2, heights=c(1,1.4)), 
          nrow=1, widths=c(2,1))

summary.plot = ggarrange(igj.plot.2a +
            scale_x_continuous(expand = expansion(mult = c(0.3, 0), 
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

set.seed(1)
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

load(paste("Results/Alveolates-", mytag, "/Alveolates.RData",sep=""))
set.3 = parallelised.runs[[1]]
load(paste("Results/Alveolates-", mytag ,"/Non-Alveolates.RData",sep=""))
set.4 = parallelised.runs[[1]]

fig.3c.1 = plotHypercube.bubbles(set.3) + theme(legend.position="none")
fig.3c.2 = plotHypercube.sampledgraph2(set.3, node.labels = FALSE,
                                       no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")
fig.3d.1 = plotHypercube.bubbles(set.4) + theme(legend.position="none")
fig.3d.2 = plotHypercube.sampledgraph2(set.4, node.labels = FALSE,
                                       no.times = TRUE, edge.label.size = 3) + theme(legend.position = "none")

fig.3x = ggarrange(ggarrange(fig.3a.1,fig.3a.2, nrow = 2, heights=c(1,2), labels = c("i","ii")) + theme(legend.position = "none"),
                  ggarrange(fig.3b.1,fig.3b.2, nrow = 2, heights=c(1,2),labels = c("i","ii")) + theme(legend.position = "none"),
                  ggarrange(fig.3c.1,fig.3c.2, nrow = 2, heights=c(1,2),labels = c("i","ii")) + theme(legend.position = "none"),
                  ggarrange(fig.3d.1,fig.3d.2, nrow = 2, heights=c(1,2),labels = c("i","ii")) + theme(legend.position = "none"),
                  nrow = 2, ncol = 2) #,labels = c("A","B","C","D"))
sf = 2
png("Results/fig-3x.png", width=500*sf, height=1000*sf, res=72*sf)
print(fig.3x)
dev.off()

fig.3x1 = ggarrange(          ggarrange(fig.3c.1,fig.3c.2, nrow = 2, heights=c(1,2),labels = c("i","ii")) + theme(legend.position = "none"),
                   ggarrange(fig.3d.1,fig.3d.2, nrow = 2, heights=c(1,2),labels = c("i","ii")) + theme(legend.position = "none"),
                   nrow = 1) #,labels = c("
sf = 2
png("Results/fig-3x1.png", width=600*sf, height=600*sf, res=72*sf)
print(fig.3x1)
dev.off()

#### Fig. 4
fig.4.dots = read.csv("FBA/MitoMammal/Results/MAX_ATP/single-double.csv")
fig.4.bars = read.csv("FBA/MitoMammal/Results/MAX_ATP/single-double-AOX.csv")

fig.4.dots = fig.4.dots[fig.4.dots$OBJ == "MAX_ATP",]
fig.4.bars = fig.4.bars[fig.4.bars$OBJ == "MAX_ATP",]

fig.4 = ggplot()+
  geom_col(data = fig.4.bars, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_normoxic), fill = "#4A4A4A", color = "white")+
  geom_point(data = fig.4.dots, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_normoxic), color = "#D9D9D9", size = 2)+
  labs(x = "Knockout", y = "Modelled ATP production")+
  theme_minimal() 

sf = 2
png("Results/fig-4.png", width = 800*sf, height = 600*sf, res = 72*sf)
fig.4 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(legend.position="none")
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


#### SI fig. 5:
fig.5.data = read.csv("FBA/MitoMammal/Results/MAX_ATP/single-double-AOX-NDH2.csv")
fig.5.atp = fig.5.data[fig.5.data$OBJ == "MAX_ATP",]
fig.5.pmf = fig.5.data[fig.5.data$OBJ == "MAX_PMF",]

fig.5.top = ggplot()+
  geom_col(data = fig.5.atp, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_normoxic), fill = "#FF000050")+
  geom_col(data = fig.5.pmf, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_normoxic/3.737), fill = "#0000FF50")+
  labs(x = "KO", y = "Normoxic max ATP or\n (scaled) PMR generation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(legend.position="none")
fig.5.bot = ggplot()+
  geom_col(data = fig.5.atp, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_hypoxic), fill = "#FF000050")+
  geom_col(data = fig.5.pmf, aes(x = factor(KO, levels = unique(KO), ordered = T), y = max_obj_hypoxic/3.737), fill = "#0000FF50")+
  labs(x = "KO", y = "Hypoxic max ATP or\n (scaled) PMR generation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(legend.position="none")

g.fig.5 = ggarrange(fig.5.top,fig.5.bot,nrow =2)

sf = 2
png("Results/SI-fig-5.png", width = 800*sf, height = 600*sf, res = 72*sf)
print(g.fig.5)
dev.off()

