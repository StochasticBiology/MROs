#!/usr/bin/env Rscript

### A script to perform systematic knockout (or knockdown) of genes in a particular GEM
## We're interested in how the knockdown/-out of a particular gene changes metabolic requirements for other genes and their products, particularly for mitochondrial functions.

### Requires the packages libSBML and sybilSBML, the former of which require installation system wide. Similarly, you need the package Rglpk, which needs glpk system wide.
## Note that there appears to be a bug for the installation on Windows for the time being (Feb 20, 2025), so I advice you use on a Linux operating system
## Note also that installing all required packages, particularly libSBML and Rglpk, requires some extra work. I recommend you consult the following sources in order to install them: https://github.com/SysBioChalmers/sybil-SBML#installation and https://stackoverflow.com/questions/32647337/installing-rglpk-with-rstudio

require(libSBML)
require(sybilSBML)
require(rsbml) 
require(Rglpk)
require(glpkAPI)
require(ggplot2)
require(ggarrange)

# Read the SBML model mitoMammal
mod <- sybilSBML::readSBMLmod("mito-mammal.xml")


# Go through each of the reaction, CI-V (NADH, SDH, COB, COX, ATP) and PDH, TCA, Fe-S cluster assembly(?), knock down/out and see what happens to the rest

## Reactions
# Objective reactions are max ATP, max heme, max lipid synthesis (IMM), and max amino acid synthesis and are numbered 71 through 74 (1-70 are exchange reactions)
OBJ_inds <- c(71:74)

# Make a separate model for each of the main objectives
mod.ATP  <- changeObjFunc(mod, react = 71) # 1 is default for obj_coef; rest is zero
mod.HEME <- changeObjFunc(mod, react = 72)
mod.LIP  <- changeObjFunc(mod, react = 73)
mod.PROT <- changeObjFunc(mod, react = 74)

# CI-V (reactions 38-42 in the spreadsheet)
ETC_inds  <- c(108:112)
ETC_names <- c("NADH dehydrogenase", "succinate dehydrogenase", "cytochrome c reductase", "cytochrome c oxidase", "ATP synthase")

### TCA (28-37)
TCA_inds  <- c(98:107)
TCA_names <- c("pyruvate dehydrogenase","citrate synthase","Aconitate hydratase","Isocitrate dehydrogenase (NAD+)","Isocitrate dehydrogenase (NADP+), mitochondrial",
	       "2-oxoglutarate dehydrogenase complex", "Succinate--CoA ligase (GDP-forming)", "Succinate--CoA ligase (ADP-forming)", "fumarate hydratase",
	       "malate dehydrogenase 2, NAD (mitochondrial)")

# Haem synthesis (384-393) and degradation (394-396)
HEME_inds  <- c(384:396)
HEME_names <- c("5'-aminolevulinate synthase 1", "Unknown", "aminolevulinate dehydratase", "hydroxymethylbilane synthase",
		"uroporphyrinogen III synthase", "uroporphyrinogen decarboxylase", "coproporphyrinogen oxidase",
	 	"protoporphyrinogen oxidase", "ferrochelatase", "Unknown", "heme oxygenase 1 or 2", "biliverdin reductase B", "biliverdin reductase A")

# KNOCKOUT: ko <- optimizeProb(model, react = i, lb = 0, ub = 0, lpdir = "max"/"min") or deleteGene(
# 50% KNOCKDOWN: kd <- optimizeProb(model, react = i, lb = ub0, ub = 500, lpdir = "max"/"min") (assuming irreversibility)

cat("Producing baseline for max ATP..\n")
# No knockdown with max ATP production
bl <- optimizeProb(mod.ATP, lpdir = "max")
bl.fd <- getFluxDist(bl)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-none.csv"
cat(paste(c("Exporting to file ",filename,"...\n\n")))
write.csv(file = paste("Results/",filename), round(bl.fd,4), row.names = F)

# Knock down CI 50% compared to it baseline activity.
ETC_ind <- ETC_inds[1]
cat("Knocking down complex I (",ETC_names[1],")...\n")
ci.kd <- optimizeProb(mod.ATP, react = ETC_ind, lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ci.fd <- getFluxDist(ci.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-CI.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename), round(ci.fd,4), row.names = F)

# Knock down CII 50% compared to its baseline activity
ETC_ind <- ETC_inds[2]
cat("Knocking down complex II (",ETC_names[2],")...\n")
cii.kd <- optimizeProb(mod.ATP, react = ETC_inds[2], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cii.fd <- getFluxDist(cii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-CII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename), round(cii.fd,4), row.names = F)

# Knock down CIII 50% compared to its baseline activity
ETC_ind <- ETC_inds[3]
cat("Knocking down complex III (",ETC_names[3],")...\n")
ciii.kd <- optimizeProb(mod.ATP, react = ETC_inds[3], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ciii.fd <- getFluxDist(ciii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-CIII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename), round(ciii.fd,4), row.names = F)

# Knock down CIV 50% compared to its baseline activity
ETC_ind <- ETC_inds[4]
cat("Knocking down complex IV (",ETC_names[4],")...\n")
civ.kd <- optimizeProb(mod.ATP, react = ETC_inds[4], lb = 0, ub = bl.fd[ETC_ind]/2, lpdir = "max")
civ.fd <- getFluxDist(civ.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-CIV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename), round(civ.fd,4), row.names = F)

# Knock down CV 50% compared to its baseline activity
ETC_ind <- ETC_inds[5]
cat("Knocking down complex V (",ETC_names[5],")...\n")
cv.kd <- optimizeProb(mod.ATP, react = ETC_inds[5], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cv.fd <- getFluxDist(cv.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-CV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename), round(cv.fd,4), row.names = F)

# Knock down PDH 50% compared to its baseline activity
cat("Knocking down PDH ...\n")
TCA_ind <- TCA_inds[1]
pdh.kd <- optimizeProb(mod.ATP, react = TCA_inds[1], lb = -bl.fd[TCA_ind]/2, ub = bl.fd[TCA_ind]/2, lpdir = "max")
pdh.fd <- getFluxDist(pdh.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "ATP-PDH.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(pdh.fd,4), row.names = F)

####################

cat("Producing baseline for max HEME..\n")
# No knockdown with max HEME production
bl <- optimizeProb(mod.HEME, lpdir = "max")
bl.fd <- getFluxDist(bl)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-none.csv"
cat(paste(c("Exporting to file ",filename,"...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(bl.fd,4), row.names = F)

# Knock down CI 50% compared to it baseline activity. Since baseline activity = 0, we use -500 and 500
ETC_ind <- ETC_inds[1]
cat("Knocking down complex I (",ETC_names[1],")...\n")
ci.kd <- optimizeProb(mod.HEME, react = ETC_inds[1], lb = -500, ub = 500, lpdir = "max")
ci.fd <- getFluxDist(ci.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-CI.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ci.fd,4), row.names = F)

# Knock down CII 50% compared to its baseline activity
ETC_ind <- ETC_inds[2]
cat("Knocking down complex II (",ETC_names[2],")...\n")
cii.kd <- optimizeProb(mod.HEME, react = ETC_inds[2], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cii.fd <- getFluxDist(cii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-CII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cii.fd,4), row.names = F)

# Knock down CIII 50% compared to its baseline activity
ETC_ind <- ETC_inds[3]
cat("Knocking down complex III (",ETC_names[3],")...\n")
ciii.kd <- optimizeProb(mod.HEME, react = ETC_inds[3], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ciii.fd <- getFluxDist(ciii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-CIII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ciii.fd,4), row.names = F)

# Knock down CIV 50% compared to its baseline activity
ETC_ind <- ETC_inds[4]
cat("Knocking down complex IV (",ETC_names[4],")...\n")
civ.kd <- optimizeProb(mod.HEME, react = ETC_inds[4], lb = 0, ub = bl.fd[ETC_ind]/2, lpdir = "max")
civ.fd <- getFluxDist(civ.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-CIV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(civ.fd,4), row.names = F)

# Knock down CV 50% compared to its baseline activity; since it is used in reverse for the base case, we reverse:
ETC_ind <- ETC_inds[5]
cat("Knocking down complex V (",ETC_names[5],")...\n")
cv.kd <- optimizeProb(mod.HEME, react = ETC_inds[5], lb = bl.fd[ETC_ind]/2, ub = -bl.fd[ETC_ind]/2, lpdir = "max")
cv.fd <- getFluxDist(cv.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-CV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cv.fd,4), row.names = F)

# Knock down PDH 50% compared to its baseline activity
cat("Knocking down PDH ...\n")
TCA_ind <- TCA_inds[1]
pdh.kd <- optimizeProb(mod.HEME, react = TCA_inds[1], lb = -bl.fd[TCA_ind]/2, ub = bl.fd[TCA_ind]/2, lpdir = "max")
pdh.fd <- getFluxDist(pdh.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "HEME-PDH.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(pdh.fd,4), row.names = F)

####################

cat("Producing baseline for max LIPID..\n")
# No knockdown with max LIPID production
bl <- optimizeProb(mod.LIP, lpdir = "max")
bl.fd <- getFluxDist(bl)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-none.csv"
cat(paste(c("Exporting to file ",filename,"...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(bl.fd,4), row.names = F)

# Knock down CI 50% compared to it baseline activity. Since baseline activity = 0, we use -500 and 500
ETC_ind <- ETC_inds[1]
cat("Knocking down complex I (",ETC_names[1],")...\n")
ci.kd <- optimizeProb(mod.LIP, react = ETC_inds[1], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ci.fd <- getFluxDist(ci.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-CI.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ci.fd,4), row.names = F)

# Knock down CII 50% compared to its baseline activity
ETC_ind <- ETC_inds[2]
cat("Knocking down complex II (",ETC_names[2],")...\n")
cii.kd <- optimizeProb(mod.LIP, react = ETC_inds[2], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cii.fd <- getFluxDist(cii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-CII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cii.fd,4), row.names = F)

# Knock down CIII 50% compared to its baseline activity
ETC_ind <- ETC_inds[3]
cat("Knocking down complex III (",ETC_names[3],")...\n")
ciii.kd <- optimizeProb(mod.LIP, react = ETC_inds[3], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ciii.fd <- getFluxDist(ciii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-CIII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ciii.fd,4), row.names = F)

# Knock down CIV 50% compared to its baseline activity
ETC_ind <- ETC_inds[4]
cat("Knocking down complex IV (",ETC_names[4],")...\n")
civ.kd <- optimizeProb(mod.LIP, react = ETC_inds[4], lb = 0, ub = bl.fd[ETC_ind]/2, lpdir = "max")
civ.fd <- getFluxDist(civ.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-CIV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(civ.fd,4), row.names = F)

# Knock down CV 50% compared to its baseline activity; since it is used in reverse for the base case, we reverse:
ETC_ind <- ETC_inds[5]
cat("Knocking down complex V (",ETC_names[5],")...\n")
cv.kd <- optimizeProb(mod.LIP, react = ETC_inds[5], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cv.fd <- getFluxDist(cv.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-CV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cv.fd,4), row.names = F)

# Knock down PDH 50% compared to its baseline activity
cat("Knocking down PDH ...\n")
TCA_ind <- TCA_inds[1]
pdh.kd <- optimizeProb(mod.LIP, react = TCA_inds[1], lb = -bl.fd[TCA_ind]/2, ub = bl.fd[TCA_ind]/2, lpdir = "max")
pdh.fd <- getFluxDist(pdh.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "LIPID-PDH.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(pdh.fd,4), row.names = F)

####################

cat("Producing baseline for max PROTEIN..\n")
# No knockdown with max PROTEIN production
bl <- optimizeProb(mod.PROT, lpdir = "max")
bl.fd <- getFluxDist(bl)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-none.csv"
cat(paste(c("Exporting to file ",filename,"...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(bl.fd,4), row.names = F)

# Knock down CI 50% compared to it baseline activity. Since baseline activity = 0, we use -500 and 500
ETC_ind <- ETC_inds[1]
cat("Knocking down complex I (",ETC_names[1],")...\n")
ci.kd <- optimizeProb(mod.PROT, react = ETC_inds[1], lb = -500, ub = 500, lpdir = "max")
ci.fd <- getFluxDist(ci.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-CI.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ci.fd,4), row.names = F)

# Knock down CII 50% compared to its baseline activity
ETC_ind <- ETC_inds[2]
cat("Knocking down complex II (",ETC_names[2],")...\n")
cii.kd <- optimizeProb(mod.PROT, react = ETC_inds[2], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
cii.fd <- getFluxDist(cii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-CII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cii.fd,4), row.names = F)

# Knock down CIII 50% compared to its baseline activity
ETC_ind <- ETC_inds[3]
cat("Knocking down complex III (",ETC_names[3],")...\n")
ciii.kd <- optimizeProb(mod.PROT, react = ETC_inds[3], lb = -bl.fd[ETC_ind]/2, ub = bl.fd[ETC_ind]/2, lpdir = "max")
ciii.fd <- getFluxDist(ciii.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-CIII.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(ciii.fd,4), row.names = F)

# Knock down CIV 50% compared to its baseline activity
ETC_ind <- ETC_inds[4]
cat("Knocking down complex IV (",ETC_names[4],")...\n")
civ.kd <- optimizeProb(mod.PROT, react = ETC_inds[4], lb = 0, ub = bl.fd[ETC_ind]/2, lpdir = "max")
civ.fd <- getFluxDist(civ.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-CIV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(civ.fd,4), row.names = F)

# Knock down CV 50% compared to its baseline activity; since it is used in reverse for the base case, we reverse:
ETC_ind <- ETC_inds[5]
cat("Knocking down complex V (",ETC_names[5],")...\n")
cv.kd <- optimizeProb(mod.PROT, react = ETC_inds[5], lb = bl.fd[ETC_ind]/2, ub = -bl.fd[ETC_ind]/2, lpdir = "max")
cv.fd <- getFluxDist(cv.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-CV.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(cv.fd,4), row.names = F)

# Knock down PDH 50% compared to its baseline activity
cat("Knocking down PDH ...\n")
TCA_ind <- TCA_inds[1]
pdh.kd <- optimizeProb(mod.PROT, react = TCA_inds[1], lb = -bl.fd[TCA_ind]/2, ub = bl.fd[TCA_ind]/2, lpdir = "max")
pdh.fd <- getFluxDist(pdh.kd)#[c(OBJ_inds,ETC_inds,TCA_inds,HEME_inds)]

filename <- "PROTEIN-PDH.csv"
cat(paste(c("Exporting to file ", filename, "...\n\n")))
write.csv(file = paste("Results/",filename, sep = ""), round(pdh.fd,4), row.names = F)