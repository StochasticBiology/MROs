#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there are exactly three arguments: if not, return an error
if (length(args)!=4) {
  stop("Need four arguments: a Phyllip or Newick tree, binary strings of loss patterns, an output folder, and an output label (1: 1b --> 0; 2: 1b --> 1).\n Usage: MRO-rcg-uncertainty.R [tree file] [barcodes] [output folder] [output label]  \n", call.=FALSE)
}

library(phytools)
library(parallel)
library(sybilSBML)
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

fba.mod = readSBMLmod("FBA/MitoMammal/MitoMammal.xml")
etc.inds <- c(98,108:112)
etc.rs = c("PDH","CI","CII","CIII","CIV","CV")#fba.mod@react_name[etc.inds]
EX_o2 <- 50
 

atp.obj = rep(0, length(fba.mod@react_name))
atp.obj[71] = 1
res.df = data.frame(name = NULL, ko.val = NULL, ko.atp.val = NULL, ko.val.hypoxia = NULL, ko.atp.val.hypoxia = NULL)
for(ko in 1:length(etc.rs)) {
  tmp.mod = fba.mod
  tmp.mod@lowbnd[etc.rs[ko]] = 0
  tmp.mod@uppbnd[etc.rs[ko]] = 0
  soln = optimizeProb(tmp.mod)
  soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
  tmp.mod@lowbnd[EX_o2] = 0
  tmp.mod@uppbnd[EX_o2] = 1
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
  res.df = rbind(res.df, data.frame(name=etc.rs[ko], 
                                    ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                    ko.val.hypoxia=soln.hypoxia@lp_obj,
                                    ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))
}

tmp.mod = fba.mod
tmp.mod@lowbnd[98:107] = 0
tmp.mod@uppbnd[98:107] = 0
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
tmp.mod@lowbnd[EX_o2] = 0
tmp.mod@uppbnd[EX_o2] = 1
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="TCA", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))

tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
tmp.mod@lowbnd[EX_o2] = 0
tmp.mod@uppbnd[EX_o2] = 1
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="full model", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))

write.csv(res.df, "fba-res.csv")
res.df = read.csv("fba-res.csv")

#print(res.df)

#ggplot(res.df, aes(x=ko.atp.val, y=ko.atp.val.hypoxia, label=label)) + geom_point() + ggrepel::geom_text_repel()

#ggplot(res.df, aes(x=ko.val, y=ko.val.hypoxia, label=label)) + 
#  geom_point() + 
#  ggrepel::geom_text_repel()
# First-loss probability vs metabolic impact from FBA in the normoxic case

# Labelling for easy access to elements
res.df$label = ""
res.df$label[1:7]  = c("PDH","CI","CII","CIII","CIV","CV","TCA")
#res.df$label[2] = "CI"
#res.df$label[3] = "CII"
#res.df$label[4] = "CIII"
#res.df$label[5] = "CIV"
#res.df$label[6] = "CV"
#res.df$label[7] = "TCA"

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
#print(first.df)
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

png(paste(outfolder, "mro-first-fba-normoxic-", tag, ".png", sep = ""), width=300*sf, height=300*sf, res=72*sf)
print(g.normoxic + ggtitle(tstr))
dev.off()

# First-loss probability vs metabolic impact from FBA in the anoxic case

# Re-labelling for easy access
res.df$label = ""
res.df$label[1:7]  = c("PDH","CI","CII","CIII","CIV","CV","TCA")
#res.df$label[2] = "CI"
#res.df$label[3] = "CII"
#res.df$label[4] = "CIII"
#res.df$label[5] = "CIV"
#res.df$label[6] = "CV"
#res.df$label[7] = "TCA"

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
#print(first.df)
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

png(paste(outfolder, "mro-first-fba-anoxic-", tag, ".png",sep = ""), width=300*sf, height=300*sf, res=72*sf)
print(g.anoxic + ggtitle(tstr))
dev.off()

save.image(paste(outfolder, "mro-save-", tag,".RData",sep = ""))