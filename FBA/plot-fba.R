#!/usr/bin/env Rscript

require(ggplot2)
require(forcats)
#require(ggarrange)

# Import data for unrestricted oxygen case (base case) for single and double KOs
d <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-no-restriction.csv")

# Import data for restricted oxygen (10% of oxygen used at at no KO) for single and double KOs
d2 <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-oxygen-restriction.csv")

# Import data for unrestricted oxygen case (base case) for single and double KOs
d3 <- read.csv("MitoMammal/Results/triple-quadruple-KO-MAX_ATP-no-restriction.csv")

# Import data for restricted oxygen (10% of oxygen used at at no KO) for single and double KOs
d4 <- read.csv("MitoMammal/Results/triple-quadruple-KO-MAX_ATP-oxygen-restriction.csv")

p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3 <- ggplot(data = d3)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p4 <- ggplot(data = d4)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
  lims(y=c(0,110))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 900
h = 500
res <- 3
filename <- "single-double-MAX_ATP-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-MAX_ATP-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

filename <- "triple-quadruple-MAX_ATP-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p3
dev.off()

filename <- "triple-quadruple-MAX_ATP-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p4
dev.off()