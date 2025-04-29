#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there are exactly 5 arguments: if not, return an error
if (length(args)!=5) {
  stop("Need 5 arguments: a Phyllip or Newick tree, binary strings of loss patterns, an output folder, the penalty to employ (0/1), whether the data is accompanied by timings (0/1), and an input label.\n Usage: runHyperTraPS-inference.R [tree file] [barcodes] [output label] [penalty (0/1)] [with_timings (0/1)]\n", call.=FALSE)
}

require(phytools)
require(parallel)
require(sybilSBML)
require(ggplot2)

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

# harmonise labels across tree and barcode data # Download the full tree without changing type of file or labels in the file itself
#tree$tip.label = gsub("'", "", tree$tip.label)
tree$tip.label = gsub("_", " ", tree$tip.label)

# TODO: A much more robust way of getting the subtrees of interest, after which curate.tree will pick out the relevant tips, outputting the message that it found observations not corresponding to tips, which are OK(?)
#if(grep(x = output, "Apicomplexans")){
#  tree <- extract.clade(tree, 96)
#}else if(grep(x = output, "Ciliophora")){
#  tree <- extract.clade(tree,86)
#}else if(grep(x = output, "Non-Apicomplexans")){
#  tree <- splitTree(tree, split = list(node = 96, bp = 1))[[1]]
#}else if(grep(x = output, "Non-Ciliophora")){
#  tree <- splitTree(tree, split = list(node = 86, bp = 1))[[1]]
#}else if(grep(x = output, "Unicellular")){
#
#}else if(grep(x = output, "HyperTraPS-CT")){


# ancestral state reconstruction and transition gathering
data.ncbi = curate.tree(tree, df, losses=TRUE)

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
  parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                                MoreArgs = list(src.data=data.ncbi,
                                                featurenames=featurenames),
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), n.seed))

}else if(with_timings == 1){

  # Run HyperTraPS-CT with the time interval [0,tau]
  # wrapper function for continuous-time HyperTraPS analysis
  parallel.ct.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs,
                      endtimes = src.data$times,
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }
  
  parallelised.runs <- mcmapply(parallel.ct.fn, seed=1:n.seed,
                                MoreArgs = list(src.data=data.ncbi,
                                                featurenames=featurenames),
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), n.seed))

  # Run HyperTraPS-CT with the time interval [0.5*tau',1.5*tau'] with tau = 0.01*tau
  # wrapper function for continuous-time HyperTraPS analysis
  parallel.ct2.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs,
                      starttimes = 0.5*0.01*src.data$times, endtimes = 1.5*0.01*src.data$times,
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }
  
  parallelised.ct.runs <- mcmapply(parallel.ct2.fn, seed=1:n.seed,
                                MoreArgs = list(src.data=data.ncbi,
                                                featurenames=featurenames),
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), n.seed))

}else{
  stop("with_timings must be 0 (no) or 1 (yes).", call. = FALSE)
}

# If timings were passed, run HyperTraPS (without timings) on the same data
if(with_timings){

  # wrapper function for HyperTraPS analysis
  parallel.fn = function(seed, src.data, featurenames) {
    return(HyperTraPS(src.data$dests, initialstates = src.data$srcs, 
                      losses = 1, length = 5, kernel = 6,
                      seed = seed, samplegap = 100, penalty = penalty,
                      featurenames = featurenames))
  }

  parallelised.discrete.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                                 MoreArgs = list(src.data=data.ncbi,
                                                 featurenames=featurenames),
                                 SIMPLIFY = FALSE,
                                 mc.cores = min(detectCores(), n.seed))
}

if(!dir.exists("Results"))dir.create("Results")
if(!dir.exists(outfolder))dir.create(outfolder)

save.image(paste(outfolder, "mro-save-", penalty, "-", with_timings, ".RData",sep = ""))