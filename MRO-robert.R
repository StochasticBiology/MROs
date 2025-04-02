#!/usr/bin/env Rscript

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
# IGJ March 2025 -- updated this file
source("hypertraps.R")

# read data
tree = read.tree("Data/list-species.nwk")
#tree = read.tree("Data/mro-ncbi-tree.phy")
df = read.csv("Data/mro-barcodes-2025.csv")
=======
treefile <- as.character(args[1])
datafile <- as.character(args[2])

tree <- read.tree(paste("Data/",treefile,sep =""))
df <- read.csv(paste("Data/",datafile, sep = ""))

# IGJ March 2025 -- updated this file
source("hypertraps.R")

#tree = read.tree("Data/mro-ncbi-tree-2025.phy")
#df = read.csv("Data/mro-barcodes-2025.csv")
>>>>>>> Stashed changes

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
<<<<<<< Updated upstream
                    seed = seed, samplegap = 10, penalty = 0.,
=======
                    seed = seed, samplegap = 100, penalty = 0.,
>>>>>>> Stashed changes
                    featurenames = featurenames))
}

# run these experiments in parallel. should take a few core minutes each
<<<<<<< Updated upstream
n.seed = 3
=======
n.seed = 5
>>>>>>> Stashed changes
featurenames = colnames(df)[2:length(colnames(df))]
parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                               MoreArgs = list(src.data=data.ncbi,
                                               featurenames=featurenames),
                               SIMPLIFY = FALSE,
                               mc.cores = min(detectCores(), n.seed))

g.data = plotHypercube.curated.tree(data.ncbi, scale.fn=geom_blank(), names=TRUE)
<<<<<<< Updated upstream
=======
g.blurb = plotHypercube.summary(parallelised.runs[[1]]) + theme(legend.position="none")
>>>>>>> Stashed changes
g.bub.1 = plotHypercube.bubbles(parallelised.runs[[1]]) + theme(legend.position="none")
g.bub.2 = plotHypercube.bubbles(parallelised.runs[[2]]) + theme(legend.position="none")
g.bub.3 = plotHypercube.bubbles(parallelised.runs[[3]]) + theme(legend.position="none")
g.sg.1 = plotHypercube.sampledgraph2(parallelised.runs[[1]], 
                                     thresh=0.05, use.arc=FALSE, 
                                     node.labels = FALSE,
                                     edge.label.size=4, no.times=TRUE) + 
  theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))

# compare outputs and check convergence
<<<<<<< Updated upstream
ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1)
# dynamics
g.sg.1
=======
sf = 2
png("Results/HyperTraPS/mro-convergence.png", width = 700*sf, height = 600*sf, res = 72*sf)
print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
dev.off()

# check summary
sf = 2
png("Results/HyperTraPS/mro-trace.png", width = 700*sf, height = 600*sf, res = 72*sf)
print(g.blurb)
dev.off()

# dynamics
#g.sg.1
>>>>>>> Stashed changes
# summary plot
g.summary = ggarrange(g.data, ggarrange(g.bub.1, g.sg.1, ncol=1, 
                                        heights=c(1,2), labels=c("B", "C")), 
                      nrow=1, widths=c(1.25,1), labels=c("A",""))
<<<<<<< Updated upstream
g.summary

sf = 2
png("mro-summary.png", width=700*sf, height=600*sf, res=72*sf)
=======
#g.summary

sf = 2
png("Results/HyperTraPS/mro-summary.png", width=700*sf, height=600*sf, res=72*sf)
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
png("mro-pathways.png", width=600*sf, height=600*sf, res=72*sf)
=======
png("Results/HyperTraPS/mro-pathways.png", width=600*sf, height=600*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.pathways)
dev.off()

fba.mod = readSBMLmod("FBA/MitoMammal/MitoMammal.xml")
etc.inds <- c(108:112)
etc.rs = fba.mod@react_name[etc.inds]
<<<<<<< Updated upstream
=======
EX_o2 <- 50
>>>>>>> Stashed changes

#cat("Reactions are...")
#print(etc.rs)
#cat("END\n\n")

atp.obj = rep(0, length(fba.mod@react_name))
atp.obj[71] = 1
res.df = data.frame(name = NULL, ko.val = NULL, ko.atp.val = NULL, ko.val.hypoxia = NULL, ko.atp.val.hypoxia = NULL)
for(ko in 1:length(fba.mod@react_name)) {
  tmp.mod = fba.mod
  tmp.mod@lowbnd[ko] = 0
  tmp.mod@uppbnd[ko] = 0
  soln = optimizeProb(tmp.mod)
  soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
<<<<<<< Updated upstream
  tmp.mod@lowbnd[50] = 0
  tmp.mod@uppbnd[50] = 1
=======
  tmp.mod@lowbnd[EX_o2] = 0
  tmp.mod@uppbnd[EX_o2] = 1
>>>>>>> Stashed changes
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
  res.df = rbind(res.df, data.frame(name=fba.mod@react_name[ko], 
                                    ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                    ko.val.hypoxia=soln.hypoxia@lp_obj,
                                    ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))
}
tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
<<<<<<< Updated upstream
tmp.mod@lowbnd[50] = 0
tmp.mod@uppbnd[50] = 1
=======
tmp.mod@lowbnd[EX_o2] = 0
tmp.mod@uppbnd[EX_o2] = 1
>>>>>>> Stashed changes
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="full model", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))

<<<<<<< Updated upstream

=======
write.csv(res.df, "fba-res.csv")
>>>>>>> Stashed changes
res.df = read.csv("fba-res.csv")
res.df$label = ""
res.df$label[108] = "CI"
res.df$label[109] = "CII"
res.df$label[110] = "CIII"
res.df$label[111] = "CIV"
res.df$label[112] = "CV"

<<<<<<< Updated upstream
write.csv(res.df, "fba-res.csv")
=======

>>>>>>> Stashed changes
ggplot(res.df, aes(x=ko.atp.val, y=ko.atp.val.hypoxia, label=label)) + geom_point() + ggrepel::geom_text_repel()

ggplot(res.df, aes(x=ko.val, y=ko.val.hypoxia, label=label)) + 
  geom_point() + 
  ggrepel::geom_text_repel()

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
print(first.df)
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

g.normoxic
fit.lm = summary(lm(propn~ko.atp.val, data=first.df))

cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

<<<<<<< Updated upstream
png("mro-first-fba.png", width=300*sf, height=300*sf, res=72*sf)
print(g.normoxic + ggtitle(tstr))
dev.off()
=======
png("Results/HyperTraPS/mro-first-fba.png", width=300*sf, height=300*sf, res=72*sf)
print(g.normoxic + ggtitle(tstr))
dev.off()

save.image("Results/HyperTraPS/mro-save.RData")
>>>>>>> Stashed changes
