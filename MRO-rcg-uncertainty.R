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

if(!dir.exists("Results"))dir.create("Results")
if(!dir.exists(outfolder))dir.create(outfolder)

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

save.image(paste(outfolder, "mro-save-", tag, "-0.RData",sep = ""))
