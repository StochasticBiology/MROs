#!/usr/bin/env Rscript

# Script to analyze HyperTraPS inferences output from runHyperTraPS-inference.R. The script accepts the argument with_timings (0/1), which, if it is zero, only produces HyperTraPS inferences without timing information. If with_timings is 1, it prints figures for a) HyperTraPS-CT inference with the interval [0,tau], b) the same inference in a classical perspective (i.e., without timing information), and ci) HyperTraPS-CT inference with the interval [0.5*tau',1.5*tau'] with tau' = 0.01*tau and cii) the same inference in a classical perspective (i.e., without timing information).

args = commandArgs(trailingOnly = T)

if (length(args)!=5) {
  stop("Need five arguments: a tree file, the corresponding file with barcodes, a folder containing the RData file with the HyperTraPS inference, its associated penalty and whether or not it accompanied timings (0/1).\n Usage: MRO-rcg.R [tree file] [RData file] [output/input folder] [penalty (0/1)] [with_timings (0/1)]\n", call.=FALSE)
}

library(phytools)
library(parallel)
library(sybilSBML)
library(ggplot2)
library(hypertrapsct)

treefile <- as.character(args[1])
datafile <- as.character(args[2])
outfolder <- as.character(args[3])
penalty <- as.numeric(args[4])
with_timings <- as.numeric(args[5])

tree <- read.tree(paste("Data/",treefile,sep =""))
df <- read.csv(paste("Data/",datafile, sep = ""))
outfolder <- paste("Results/",outfolder,"/",sep = "")
if(!dir.exists(outfolder))stop("Need an input folder with an RData file containing a HyperTraPS inference!", call. = FALSE)
load(paste(outfolder, "mro-save-", penalty, "-", with_timings, ".RData", sep = ""))

#print(outfolder)
# Harmonize labels
tree$tip.label = gsub("_", " ", tree$tip.label)

# ancestral state reconstruction and transition gathering
data.ncbi = curate.tree(tree, df, losses=TRUE)

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
png(paste(outfolder, "mro-convergence-", penalty, "-", with_timings, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
dev.off()

# check trace
sf = 2
png(paste(outfolder, "mro-trace-", penalty, "-", with_timings, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
print(g.blurb)
dev.off()

# check summary
g.summary = ggarrange(g.data, ggarrange(g.bub.1, g.sg.1, ncol=1, 
                                        heights=c(1,2), labels=c("B", "C")), 
                      nrow=1, widths=c(1.25,1), labels=c("A",""))

sf = 2
png(paste(outfolder, "mro-summary-", penalty, "-", with_timings, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
print(g.summary)
dev.off()

# Check influences
g.influences <- plotHypercube.influences(parallelised.runs[[1]])
g.influencegraph <- plotHypercube.influencegraph(parallelised.runs[[1]])

g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))

sf = 2
png(paste(outfolder, "mro-influences-", penalty, "-", with_timings, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
print(g.infl)
dev.off()

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

sf = 2
png(paste(outfolder, "mro-pathways-", penalty, "-", with_timings, ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
print(g.pathways)
dev.off()

fba.mod = readSBMLmod("FBA/MitoMammal/MitoMammal.xml")
etc.inds <- c(98,108:112)
etc.rs = c("PDH","CI","CII","CIII","CIV","CV")#fba.mod@react_name[etc.inds]
EX_o2 <- 50
TCA_inds <- c(98:107)

# Add AOX with lb = ub = 0
fba.mod = addReact(fba.mod, id = "AOX", met = c("q10h2[m]","o2[m]","q10[m]","h2o[m]"), Scoef = c(-2,-1,2,2), reversible = F, lb = 0, ub = 0)
aox.ind = 561
aox.lowbnd = 0
aox.uppbnd = 1000

atp.obj = rep(0, length(fba.mod@react_name))
atp.obj[71] = 1
res.df = data.frame(name = NULL, ko.val = NULL, ko.atp.val = NULL, ko.val.hypoxia = NULL, ko.atp.val.hypoxia = NULL)

tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
f.o2.l = 0# = -getFluxDist(soln)[EX_o2]
f.o2.u = 1
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="full model", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))

for(ko in 1:length(etc.rs)) {
  tmp.mod = fba.mod
  tmp.mod@lowbnd[etc.inds[ko]] = 0
  tmp.mod@uppbnd[etc.inds[ko]] = 0
  KOname = etc.rs[ko]
  
  # If CI is KOed, add NDH2
  #if(KOname == "CI"){
  #  tmp.mod = addReact(tmp.mod, id = "NDH2", met = c("h[m]","nadh[m]","q10[m]","o2[m]","nad[m]","q10h2[m]","o2s[m]"),
  #                     Scoef = c(-1,-1,-.999,-.002,1,.999,.002), reversible = F, lb = 0, ub = 1000)
  #}  
  if(KOname == "CIII" | KOname == "CIV"){
    tmp.mod@lowbnd[aox.ind] = aox.lowbnd
    tmp.mod@uppbnd[aox.ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
  res.df = rbind(res.df, data.frame(name=etc.rs[ko], 
                                    ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                    ko.val.hypoxia=soln.hypoxia@lp_obj,
                                    ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))
}

tmp.mod = fba.mod
tmp.mod@lowbnd[TCA_inds] = 0
tmp.mod@uppbnd[TCA_inds] = 0
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="TCA", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))

write.csv(res.df, "FBA/fba-res.csv")
res.df = read.csv("FBA/fba-res.csv")

# First-loss probability vs metabolic impact from FBA in the normoxic case

# Labelling for easy access to elements
res.df$label = ""
res.df$label[1:7]  = c("PDH","CI","CII","CIII","CIV","CV","TCA")

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
first.df$propn = 0
for(i in 1:nrow(first.df)) {
  ref = which(featurenames == first.df$label[i])-1
  propn = sum(first.steps == ref) / length(first.steps)
  first.df$propn[i] = propn
}
g.normoxic = ggplot(first.df, aes(x=ko.atp.val, y=propn, label=label)) + 
  geom_point() + geom_smooth(method="lm", color="#AAAAFF", fill="#AAAAFF") + 
  ggrepel::geom_text_repel() + theme_light() +
  labs(x = "ATP objective on KO", y = "First loss probability")

#g.normoxic
fit.lm = summary(lm(propn~ko.atp.val, data=first.df))
cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

png(paste(outfolder, "mro-first-fba-normoxic-", penalty, "-", with_timings, ".png", sep = ""), width=300*sf, height=300*sf, res=72*sf)
print(g.normoxic + ggtitle(tstr))
dev.off()

# First-loss probability vs metabolic impact from FBA in the anoxic case

# Re-labelling for easy access
res.df$label = ""
res.df$label[1:7]  = c("PDH","CI","CII","CIII","CIV","CV","TCA")

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
first.df$propn = 0
for(i in 1:nrow(first.df)) {
  ref = which(featurenames == first.df$label[i])-1
  propn = sum(first.steps == ref) / length(first.steps)
  first.df$propn[i] = propn
}
g.anoxic = ggplot(first.df, aes(x=ko.atp.val.hypoxia, y=propn, label=label)) + 
  geom_point() + geom_smooth(method="lm", color="#AAAAFF", fill="#AAAAFF") + 
  ggrepel::geom_text_repel() + theme_light() +
  labs(x = "ATP objective on KO", y = "First loss probability")

#g.anoxic
fit.lm = summary(lm(propn~ko.atp.val.hypoxia, data=first.df))
cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

png(paste(outfolder, "mro-first-fba-anoxic-", penalty, "-", with_timings, ".png",sep = ""), width=300*sf, height=300*sf, res=72*sf)
print(g.anoxic + ggtitle(tstr))
dev.off()

# If timings were passed, analyze the same data in a classical HyperTraPS picture
if(with_timings){
  
  load(paste(outfolder, "mro-save-0-1.RData", sep = ""))
  
  g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
  g.blurb = plotHypercube.summary(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.1 = plotHypercube.bubbles(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.2 = plotHypercube.bubbles(parallelised.discrete.runs[[2]]) + theme(legend.position="none")
  g.bub.3 = plotHypercube.bubbles(parallelised.discrete.runs[[3]]) + theme(legend.position="none")
  g.sg.1 = plotHypercube.sampledgraph2(parallelised.discrete.runs[[1]], 
                                       thresh=0.05, use.arc=FALSE, 
                                       node.labels = FALSE,
                                       edge.label.size=4, no.times=TRUE) + 
    theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))
  
  # compare outputs and check convergence
  sf = 2
  png(paste(outfolder,"mro-discrete-convergence-", penalty, "-", with_timings, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
  print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
  dev.off()
  
  # check summary
  sf = 2
  png(paste(outfolder,"mro-discrete-trace-", penalty, "-", with_timings, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
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
  png(paste(outfolder, "mro-discrete-summary-", penalty, "-", with_timings, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.summary)
  dev.off()
  
  # Pathways
  r = parallelised.discrete.runs[[1]]$routes
  r.c34tca = r[r[,1] == 2 | r[,1] == 3 | r[,1] == 7,]
  r.c1 = r[r[,1]==0,]
  
  tmp.c1 = tmp.c34tca = parallelised.runs[[1]]
  tmp.c34tca$routes = r.c34tca
  tmp.c1$routes = r.c1
  g.pathways = ggarrange(plotHypercube.motifs(parallelised.discrete.runs[[1]])+theme(legend.position="none"),
                         plotHypercube.motifs(tmp.c1)+theme(legend.position="none"),
                         plotHypercube.motifs(tmp.c34tca)+theme(legend.position="none"),
                         ncol=1, labels=c("i", "ii", "iii"))
  
  # r = parallelised.discrete.runs[[1]]$routes
  # r.ci = r[r[,1]==0,]
  # r.tca = r[r[,1]==7,]
  # r.c34 = r[r[,1]==2 | r[,1]==3,]
  # r.other = r[r[,1]!=0,]
  # r.pdh = r[r[,1] == 5,]
  # 
  # tmp.ci = tmp.other = tmp.tca = tmp.c34 = tmp.pdh = parallelised.discrete.runs[[1]]
  # tmp.ci$routes = r.ci
  # tmp.tca$routes = r.tca
  # tmp.c34$routes = r.c34
  # tmp.other$routes = r.other
  # tmp.pdh$routes = r.pdh
  # g.pathways = ggarrange(plotHypercube.motifs(parallelised.discrete.runs[[1]])+theme(legend.position="none"),
  #                        plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
  #                        plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
  #                        plotHypercube.motifs(tmp.c34)+theme(legend.position="none"),
  #                        plotHypercube.motifs(tmp.pdh)+theme(legend.position="none"), 
  #                        ncol=1, labels=c("A", "B", "C", "D", "E"))
  
  sf = 2
  png(paste(outfolder, "mro-discrete-pathways-", penalty, "-", with_timings, ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
  print(g.pathways)
  dev.off()
  
  # Check influences
  g.influences <- plotHypercube.influences(parallelised.discrete.runs[[1]])
  g.influencegraph <- plotHypercube.influencegraph(parallelised.discrete.runs[[1]])
  
  g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))
  
  sf = 2
  png(paste(outfolder, "mro-discrete-influences-", penalty, "-", with_timings, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.infl)
  dev.off()
  
  # Repeat for the alternative timing interval [0.5*tau',1.5*tau'], where tau' = 0.01*tau
  load(paste(outfolder, "mro-save-0-1-0.5-1.5.RData", sep = ""))
  
  g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
  g.blurb = plotHypercube.summary(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.1 = plotHypercube.bubbles(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.2 = plotHypercube.bubbles(parallelised.discrete.runs[[2]]) + theme(legend.position="none")
  g.bub.3 = plotHypercube.bubbles(parallelised.discrete.runs[[3]]) + theme(legend.position="none")
  g.sg.1 = plotHypercube.sampledgraph2(parallelised.discrete.runs[[1]], 
                                       thresh=0.05, use.arc=FALSE, 
                                       node.labels = FALSE,
                                       edge.label.size=4, no.times=TRUE) + 
    theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))
  
  # compare outputs and check convergence
  sf = 2
  png(paste(outfolder,"mro-convergence-", penalty,"-", with_timings, "-0.5-1.5.png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
  print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
  dev.off()
  
  # check summary
  sf = 2
  png(paste(outfolder,"mro-trace-", penalty,"-", with_timings, "-0.5-1.5.png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
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
  png(paste(outfolder, "mro-summary-", penalty,"-", with_timings, "-0.5-1.5.png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.summary)
  dev.off()
  
  # Pathways
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
  
  sf = 2
  png(paste(outfolder, "mro-pathways-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
  print(g.pathways)
  dev.off()
  
  # Check influences
  g.influences <- plotHypercube.influences(parallelised.discrete.runs[[1]])
  g.influencegraph <- plotHypercube.influencegraph(parallelised.discrete.runs[[1]])
  
  g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))
  
  sf = 2
  png(paste(outfolder, "mro-influences-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.infl)
  dev.off()
  
  load(paste(outfolder, "mro-save-0-1-0.5-1.5.RData", sep = ""))
  
  g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
  g.blurb = plotHypercube.summary(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.1 = plotHypercube.bubbles(parallelised.discrete.runs[[1]]) + theme(legend.position="none")
  g.bub.2 = plotHypercube.bubbles(parallelised.discrete.runs[[2]]) + theme(legend.position="none")
  g.bub.3 = plotHypercube.bubbles(parallelised.discrete.runs[[3]]) + theme(legend.position="none")
  g.sg.1 = plotHypercube.sampledgraph2(parallelised.discrete.runs[[1]], 
                                       thresh=0.05, use.arc=FALSE, 
                                       node.labels = FALSE,
                                       edge.label.size=4, no.times=TRUE) + 
    theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))
  
  # compare outputs and check convergence
  sf = 2
  png(paste(outfolder,"mro-discrete-convergence-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
  print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
  dev.off()
  
  # check summary
  sf = 2
  png(paste(outfolder,"mro-discrete-trace-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
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
  png(paste(outfolder, "mro-discrete-summary-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.summary)
  dev.off()
  
  # Pathways
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
  
  sf = 2
  png(paste(outfolder, "mro-discrete-pathways-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
  print(g.pathways)
  dev.off()
  
  # Check influences
  g.influences <- plotHypercube.influences(parallelised.discrete.runs[[1]])
  g.influencegraph <- plotHypercube.influencegraph(parallelised.discrete.runs[[1]])
  
  g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))
  
  sf = 2
  png(paste(outfolder, "mro-discrete-influences-", penalty, "-", with_timings, "-0.5-1.5.png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.infl)
  dev.off()
  
}