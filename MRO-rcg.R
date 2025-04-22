#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there are exactly five arguments: if not, return an error
if (length(args)!=5) {
<<<<<<< Updated upstream
  stop("Need five arguments: a Phyllip or Newick tree, binary strings of loss patterns, an output folder, the penalty to employ (0/1), an whether the data is accompanied by timings (0/1).\n Usage: MRO-rcg.R [tree file] [barcodes] [output label] [penalty (0/1)] [with_timings (0/1)] \n", call.=FALSE)
=======
  stop("Need five arguments: a Phyllip or Newick tree, binary strings of loss patterns, an output folder, the penalty to employ (0/1), whether the data is accompanied by timings (0/1).\n Usage: MRO-rcg.R [tree file] [barcodes] [output label] [penalty (0/1)] [with_timings (0/1)]\n", call.=FALSE)
>>>>>>> Stashed changes
}

library(phytools)
library(parallel)
library(sybilSBML)
library(ggplot2)

treefile <- as.character(args[1])
datafile <- as.character(args[2])
output <- as.character(args[3])
penalty <- as.numeric(args[4])
with_timings <- as.numeric(args[5])

tree <- read.tree(paste("Data/",treefile,sep =""))
df <- read.csv(paste("Data/",datafile, sep = ""))
outfolder <- paste("Results/",output,"/", sep = "")

# IGJ March 2025 -- updated this file
source("hypertraps.R")

# harmonise labels across tree and barcode data
# IGJ March 2025 -- think I added something here
# df$Organism = gsub("_", " ", df$Organism)

#tree$tip.label = gsub("'", "", tree$tip.label)
tree$tip.label = gsub("_", " ", tree$tip.label)

# ancestral state reconstruction and transition gathering
data.ncbi = curate.tree(tree, df, losses=TRUE)

<<<<<<< Updated upstream
# wrapper function for HyperTraPS analysis
parallel.fn = function(seed, src.data, featurenames) {
  return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                    losses = 1, length = 5, kernel = 6,
                    seed = seed, samplegap = 100, penalty = penalty,
                    featurenames = featurenames))
}

# wrapper function for continuous-time HyperTraPS analysis
parallel.ct.fn = function(seed, src.data, featurenames) {
  return(HyperTraPS(src.data$dests, initialstates = src.data$srcs,
		    endtimes = src.data$times, 
                    losses = 1, length = 5, kernel = 6,
                    seed = seed, samplegap = 100, penalty = penalty,
                    featurenames = featurenames))
}

# run these experiments in parallel. should take a few core minutes each
n.seed = 5
featurenames = colnames(df)[2:length(colnames(df))]
if(with_timings == 0){
=======
print(data.ncbi$data)

# run experiments in parallel. should take a few core minutes each
n.seed = 5
featurenames = colnames(df)[2:length(colnames(df))]

if(with_timings == 0){
  # wrapper function for HyperTraPS analysis
  parallel.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }
>>>>>>> Stashed changes
  parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                                MoreArgs = list(src.data=data.ncbi,
                                                featurenames=featurenames),
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), n.seed))
<<<<<<< Updated upstream
}else if(with_timings == 1){
  parallelised.ct.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                                   MoreArgs = list(src.data=data.ncbi,
                                                   featurenames=featurenames),
                                   SIMPLIFY = FALSE,
                                   mc.cores = min(detectCores(), n.seed))
  parallelised.runs <- parallelised.ct.runs
=======

}else if(with_timings == 1){

  # Run HyperTraPS-CT with

  # wrapper function for continuous-time HyperTraPS analysis
  parallel.ct.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs,
  	  	      endtimes = src.data$times, starttimes = src.data$times,
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }
  
  parallelised.runs <- mcmapply(parallel.ct.fn, seed=1:n.seed,
                                MoreArgs = list(src.data=data.ncbi,
                                                featurenames=featurenames),
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), n.seed))

>>>>>>> Stashed changes
}else{
  stop("with_timings must be 0 (no) or 1 (yes).", call. = FALSE)
}

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
<<<<<<< Updated upstream
png(paste(outfolder, "mro-convergence-", penalty, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
=======
png(paste(outfolder, "mro-convergence-", penalty, "-instant", ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
>>>>>>> Stashed changes
print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
dev.off()

# check summary
sf = 2
<<<<<<< Updated upstream
png(paste(outfolder, "mro-trace-", penalty, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
=======
png(paste(outfolder, "mro-trace-", penalty, "-instant", ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
png(paste(outfolder, "mro-summary-", penalty, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
=======
png(paste(outfolder, "mro-summary-", penalty, "-instant", ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.summary)
dev.off()

# Check influences
g.influences <- plotHypercube.influences(parallelised.runs[[1]])
g.influencegraph <- plotHypercube.influencegraph(parallelised.runs[[1]])

g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))

sf = 2
<<<<<<< Updated upstream
png(paste(outfolder, "mro-influences-", penalty, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
=======
png(paste(outfolder, "mro-influences-", penalty, "-instant", ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.infl)
dev.off()

r = parallelised.runs[[1]]$routes
r.ci = r[r[,1]==0,]
r.tca = r[r[,1]==7,]
r.c34 = r[r[,1]==2 | r[,1]==3,]
r.pdh = r[r[,1]==5,]
r.other = r[r[,1]!=0,]

tmp.ci = tmp.other = tmp.tca = tmp.c34 = tmp.pdh = parallelised.runs[[1]]
tmp.ci$routes = r.ci
tmp.tca$routes = r.tca
tmp.c34$routes = r.c34
tmp.other$routes = r.other
tmp.pdh$routes = r.pdh
g.pathways = ggarrange(plotHypercube.motifs(parallelised.runs[[1]])+theme(legend.position="none"),
          plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.c34)+theme(legend.position="none"),
          plotHypercube.motifs(tmp.pdh)+theme(legend.position="none"),
          ncol=1, labels=c("A", "B", "C", "D", "E"))

sf = 2
<<<<<<< Updated upstream
png(paste(outfolder, "mro-pathways-", penalty, ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
=======
png(paste(outfolder, "mro-pathways-", penalty, "-instant", ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.pathways)
dev.off()

fba.mod = readSBMLmod("FBA/MitoMammal/MitoMammal.xml")
etc.inds <- c(98,108:112)
etc.rs = fba.mod@react_name[etc.inds]
EX_o2 <- 50

<<<<<<< Updated upstream
#cat("Reactions are...")
#print(etc.rs)
#cat("END\n\n")
=======
>>>>>>> Stashed changes

atp.obj = rep(0, length(fba.mod@react_name))
atp.obj[71] = 1
res.df = data.frame(name = NULL, ko.val = NULL, ko.atp.val = NULL, ko.val.hypoxia = NULL, ko.atp.val.hypoxia = NULL)
for(ko in 1:length(fba.mod@react_name)) {
  tmp.mod = fba.mod
  tmp.mod@lowbnd[ko] = 0
  tmp.mod@uppbnd[ko] = 0
  soln = optimizeProb(tmp.mod)
  soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
  tmp.mod@lowbnd[EX_o2] = 0
  tmp.mod@uppbnd[EX_o2] = 1
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
  res.df = rbind(res.df, data.frame(name=fba.mod@react_name[ko], 
                                    ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                    ko.val.hypoxia=soln.hypoxia@lp_obj,
                                    ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))
}

tmp.mod = fba.mod
tmp.mod@lowbnd[99:108] = 0
tmp.mod@uppbnd[99:108] = 0
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
res.df$label = ""
res.df$label[98]  = "PDH"
#res.df$label[99]  = "CS"
#res.df$label[100]  = "ACONT"
#res.df$label[101]  = "ICDHx"
#res.df$label[102]  = "ICDHy"
#res.df$label[103]  = "AKGD"
#res.df$label[104]  = "SUCOAS1"
#res.df$label[105]  = "SUCOAS"
#res.df$label[106]  = "FUM"
#res.df$label[107]  = "MDH"
res.df$label[108] = "CI"
res.df$label[109] = "CII"
res.df$label[110] = "CIII"
res.df$label[111] = "CIV"
res.df$label[112] = "CV"
res.df$label[561] = "TCA"

#print(res.df)

<<<<<<< Updated upstream
ggplot(res.df, aes(x=ko.atp.val, y=ko.atp.val.hypoxia, label=label)) + geom_point() + ggrepel::geom_text_repel()

ggplot(res.df, aes(x=ko.val, y=ko.val.hypoxia, label=label)) + 
  geom_point() + 
  ggrepel::geom_text_repel()
=======
#ggplot(res.df, aes(x=ko.atp.val, y=ko.atp.val.hypoxia, label=label)) + geom_point() + ggrepel::geom_text_repel()

#ggplot(res.df, aes(x=ko.val, y=ko.val.hypoxia, label=label)) + 
#  geom_point() + 
#  ggrepel::geom_text_repel()
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
g.normoxic
=======
#g.normoxic
>>>>>>> Stashed changes
fit.lm = summary(lm(propn~ko.atp.val, data=first.df))

cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

<<<<<<< Updated upstream
png(paste(outfolder, "mro-first-fba-normoxic-", penalty, ".png", sep = ""), width=300*sf, height=300*sf, res=72*sf)
=======
png(paste(outfolder, "mro-first-fba-normoxic-", penalty, "-instant", ".png", sep = ""), width=300*sf, height=300*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.normoxic + ggtitle(tstr))
dev.off()

# First-loss probability vs metabolic impact from FBA in the anoxic case
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

<<<<<<< Updated upstream
g.anoxic
=======
#g.anoxic
>>>>>>> Stashed changes
fit.lm = summary(lm(propn~ko.atp.val.hypoxia, data=first.df))

cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

<<<<<<< Updated upstream
png(paste(outfolder, "mro-first-fba-anoxic-", penalty, ".png",sep = ""), width=300*sf, height=300*sf, res=72*sf)
=======
png(paste(outfolder, "mro-first-fba-anoxic-", penalty, "-instant", ".png",sep = ""), width=300*sf, height=300*sf, res=72*sf)
>>>>>>> Stashed changes
print(g.anoxic + ggtitle(tstr))
dev.off()

# If timings were passed, run HyperTraPS (without timings) on the same data
if(with_timings){
<<<<<<< Updated upstream
  # run these experiments in parallel. should take a few core minutes each
  n.seed = 5
  featurenames = colnames(df)[2:length(colnames(df))]
=======

  # wrapper function for HyperTraPS analysis
  parallel.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }

>>>>>>> Stashed changes
  parallelised.discrete.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                                 MoreArgs = list(src.data=data.ncbi,
                                                 featurenames=featurenames),
                                 SIMPLIFY = FALSE,
                                 mc.cores = min(detectCores(), n.seed))
  
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
  png(paste(outfolder,"mro-discrete-convergence-", penalty, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
  print(ggarrange(g.data, ggarrange(g.bub.1,g.bub.2,g.bub.3,nrow=3), nrow=1))
  dev.off()
  
  # check summary
  sf = 2
  png(paste(outfolder,"mro-discrete-trace-", penalty, ".png", sep = ""), width = 700*sf, height = 600*sf, res = 72*sf)
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
  png(paste(outfolder, "mro-discrete-summary-", penalty, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.summary)
  dev.off()
  
  r = parallelised.discrete.runs[[1]]$routes
  r.ci = r[r[,1]==0,]
  r.tca = r[r[,1]==7,]
  r.c34 = r[r[,1]==2 | r[,1]==3,]
  r.other = r[r[,1]!=0,]
  r.pdh = r[r[,1] == 5,]
  
  tmp.ci = tmp.other = tmp.tca = tmp.c34 = tmp.pdh = parallelised.discrete.runs[[1]]
  tmp.ci$routes = r.ci
  tmp.tca$routes = r.tca
  tmp.c34$routes = r.c34
  tmp.other$routes = r.other
  tmp.pdh$routes = r.pdh
  g.pathways = ggarrange(plotHypercube.motifs(parallelised.discrete.runs[[1]])+theme(legend.position="none"),
            plotHypercube.motifs(tmp.ci)+theme(legend.position="none"),
            plotHypercube.motifs(tmp.tca)+theme(legend.position="none"),
            plotHypercube.motifs(tmp.c34)+theme(legend.position="none"),
            plotHypercube.motifs(tmp.pdh)+theme(legend.position="none"), 
            ncol=1, labels=c("A", "B", "C", "D", "E"))
  
  sf = 2
  png(paste(outfolder, "mro-discrete-pathways-", penalty, ".png", sep = ""), width=600*sf, height=600*sf, res=72*sf)
  print(g.pathways)
  dev.off()
  
  # Check influences
  g.influences <- plotHypercube.influences(parallelised.discrete.runs[[1]])
  g.influencegraph <- plotHypercube.influencegraph(parallelised.discrete.runs[[1]])
  
  g.infl <- ggarrange(g.influences, g.influencegraph, nrow = 1, labels=c("A","B"))
  
  sf = 2
  png(paste(outfolder, "mro-discrete-influences-", penalty, ".png", sep = ""), width=700*sf, height=600*sf, res=72*sf)
  print(g.infl)
  dev.off()
}

save.image(paste(outfolder, "mro-save-", penalty, ".RData",sep = ""))