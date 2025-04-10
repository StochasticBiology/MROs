#!/usr/bin/env Rscript

require(ggplot2)
require(forcats)
#require(ggarrange)

# Import data for unrestricted oxygen case (base case) for single and double KOs
<<<<<<< Updated upstream
d <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-no-restriction.csv")

# Import data for restricted oxygen (10% of oxygen used at at no KO) for single and double KOs
d2 <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-oxygen-restriction.csv")

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
## Import data for unrestricted oxygen case (base case) for triple and quadruple KOs
#d3 <- read.csv("MitoMammal/Results/triple-quadruple-KO-MAX_ATP-no-restriction.csv")
#
## Import data for restricted oxygen (10% of oxygen used at at no KO) for triple and quadruple KOs
#d4 <- read.csv("MitoMammal/Results/triple-quadruple-KO-MAX_ATP-oxygen-restriction.csv")
=======
d <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-normoxic.csv")

# Import data for restricted oxygen (10% of oxygen used at at no KO) for single and double KOs
d2 <- read.csv("MitoMammal/Results/single-double-KO-MAX_ATP-anoxic.csv")
>>>>>>> Stashed changes


# Plot MAX_ATP against KOs, normalized to 100% with no KO and normoxic conditions
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP/MAX_ATP[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP/d$MAX_ATP[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

<<<<<<< Updated upstream
#p3 <- ggplot(data = d3)+
#  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
#  lims(y=c(0,110))+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
#p4 <- ggplot(data = d4)+
#  geom_col(aes(x = fct_inorder(KO), y = MAX_ATP))+
#  lims(y=c(0,110))+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
res <- 3
filename <- "single-double-MAX_ATP-no-restriction.png"
=======

w <- 1000
h = 800
res <- 3
filename <- "single-double-MAX_ATP-normoxic.png"
>>>>>>> Stashed changes
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

<<<<<<< Updated upstream
filename <- "single-double-MAX_ATP-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
filename <- "triple-quadruple-MAX_ATP-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p3
dev.off()

filename <- "triple-quadruple-MAX_ATP-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p4
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
#filename <- "triple-quadruple-MAX_ATP-no-restriction.png"
#png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
#p3
#dev.off()
#
#filename <- "triple-quadruple-MAX_ATP-oxygen-restriction.png"
#png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
#p4
#dev.off()

# Plot complexes' influence on each other (i.e., how does KOing CI influence CII, etc.)

# Plot CI against KOs for MAX_ATP
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = CI/CI[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = CI/d$CI[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
res <- 3
filename <- "single-double-CI-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-CI-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

# Plot CII against KOs for MAX_ATP
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = CII/CII[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = CII/d$CII[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
res <- 3
filename <- "single-double-CII-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-CII-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

# Plot CIII against KOs for MAX_ATP
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = CIII/CIII[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = CIII/d$CIII[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
res <- 3
filename <- "single-double-CIII-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-CIII-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

# Plot CIV against KOs for MAX_ATP
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = CIV/CIV[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = CIV/d$CIV[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
res <- 3
filename <- "single-double-CIV-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-CIV-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

# Plot CV against KOs for MAX_ATP
p <- ggplot(data = d)+
  geom_col(aes(x = fct_inorder(KO), y = CV/CV[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(data = d2)+
  geom_col(aes(x = fct_inorder(KO), y = CV/d$CV[1]))+
  lims(y=c(-.1,1.1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

w <- 700
h = 600
res <- 3
filename <- "single-double-CV-no-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

filename <- "single-double-CV-oxygen-restriction.png"
=======
filename <- "single-double-MAX_ATP-anoxic.png"
>>>>>>> Stashed changes
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
dev.off()

# Make a heatmap showing the complexes, with (NONE, KOs) on the y-axis and values of (NONE, KOs) after KO for MAX_ATP in normoxic and anoxic conditions
d.norm <- d[1:7,3:9]
d.anox <- d2[1:7,3:9]



for(i in 1:ncol(d.norm)){
<<<<<<< Updated upstream
d.anox[,i]<-as.numeric(round(d.anox[,i]/d.norm[1,i]*100,2))
d.norm[,i]<-as.numeric(round(d.norm[,i]/d.norm[1,i]*100,2))
=======
d.anox[,i]<-as.numeric(round(d.anox[,i]/d.norm[1,i],2))
d.norm[,i]<-as.numeric(round(d.norm[,i]/d.norm[1,i],2))
>>>>>>> Stashed changes
}

row.names(d.norm) <- c("NONE","CI", "CII", "CIII", "CIV", "CV", "PDH")
row.names(d.anox) <- row.names(d.norm)

grid <- expand.grid(Y=names(d.norm),X=row.names(d.norm))
for(i in 1:length(grid$X)){
  grid$Z[i] = d.norm[grid$X[i],grid$Y[i]]
}

p <- ggplot(data = grid, aes(x = X, y = Y, fill = Z))+
  geom_tile() +
<<<<<<< Updated upstream
  lims(fill = c(-10,110)) +
=======
  lims(fill = c(-1,1.1)) +
>>>>>>> Stashed changes
  scale_color_gradient2()

for(i in 1:length(grid$X)){
  grid$Z[i] = d.anox[grid$X[i],grid$Y[i]]
}

p2 <- ggplot(data = grid, aes(x = X, y = Y, fill = Z))+
  geom_tile() +
<<<<<<< Updated upstream
  lims(fill = c(-10,110)) +
=======
  lims(fill = c(-1,1.1)) +
>>>>>>> Stashed changes
  scale_color_gradient2()

w <- 700
h = 600
res <- 3
<<<<<<< Updated upstream
filename <- "KO-heatmap-no-restriction.png"
=======
filename <- "KO-heatmap-normoxic.png"
>>>>>>> Stashed changes
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p
dev.off()

<<<<<<< Updated upstream
filename <- "KO-heatmap-oxygen-restriction.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
filename <- "KO-heatmap-anoxic.png"
png(file = paste("MitoMammal/Results/",filename,sep=""), width = w*res, height = h*res, res = 72*res)
p2
>>>>>>> Stashed changes
dev.off()