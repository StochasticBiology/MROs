#!/usr/bin/env Rscript

# A short script to make n datasets in which 10% of observations in the original dataset have one entry changed from zero to one, or vice versa

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  stop("Need three arguments: The number of datasets to make, the original dataset, and an output label (1: 1b --> 0; 2: 1b --> 1)].\n Usage: randomizeData.R [integer] [dataset] [label], e.g., randomizeData.R 10 mro-barcodes-2025.csv\n")
}

n <- as.numeric(args[1])
datafile <- as.character(args[2])
tag <- as.numeric(args[3])

df <- read.csv(paste("Data/",datafile, sep = ""))
if(!dir.exists(paste("Data/RandomizedData", tag, sep = "-")))dir.create(paste("Data/RandomizedData", tag, sep = "-"))

# Set a seed
set.seed(1234)

for(k in 1:n){
  this_df <- df

  # Change an observation for 10% of the observations to the opposite (0s to 1s and 1s to 0s)
  change <- sample(1:nrow(this_df), size = round(nrow(this_df)/10))

  tmp_df <- this_df[change,]
  for(i in 1:nrow(tmp_df)){
    vec <- tmp_df[i,]
    elem <- sample(2:length(vec), size = 1)
    if(vec[elem] == 0)vec[elem] = 1
    else vec[elem]= 0
    tmp_df[i,] <- vec
  }

  this_df[change,] = tmp_df

  # Save this_df
  write.csv(this_df, file = paste("Data/RandomizedData-",tag, "/randomized-dataset-", k, ".csv",sep = ""), quote = F, row.names = F)
}