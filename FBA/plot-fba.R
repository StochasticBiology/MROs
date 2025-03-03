#!/usr/bin/env Rscript

require(ggplot2)
require(forcats)
#require(ggarrange)

# Import data for unrestricted oxygen case (base case)
d <- read.csv("MitoMammal/Results/MAX_ATP-no-restriction.csv")

# Import data for restricted oxygen (10% of oxygen used at at no KO)
d2 <- read.csv("MitoMammal/Results/MAX_ATP-oxygen-restriction.csv")

p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))

res <- 3
filename <- "MAX_ATP-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = 900*res, height = 500*res, res = 72*res)
p
dev.off()

filename <- "MAX_ATP-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = 900*res, height = 500*res, res = 72*res)
p2
dev.off()