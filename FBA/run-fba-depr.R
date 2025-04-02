#!/usr/bin/env Rscript

# This is the working directory
require(libSBML)
require(sybilSBML) # Loads sybil
require(rsbml)
require(glpkAPI)
require(Rglpk)
require(combinat)

# Load up the MitoMammal GEM
mm <- readSBMLmod("MitoMammal/MitoMammal.xml")

# Read reactions
reacts <- mm.react_name #read.csv("MitoMammal/MitoMammal-reactions.txt")

names <- reacts[,2]

# Objective indices
OBJ_inds <- c(71:74)
OBJ_names <- names[OBJ_inds]#c("OF_ATP_MitoCore","OF_HEME_MitoCore","OF_LIPID_MitoCore","OF_PROTEIN_MitoCore")

# Oxygen exchange:
EX_o2 <- 50

# Test that it's really oxygen (check)
# changeObjFunc(mm, react = EX_o2)

# Tricarboxylic acid cycle complexes
TCA_inds <- c(98:107)
TCA_names <- names[TCA_inds]#c("pyruvate dehydrogenase","citrate synthase","Aconitate hydratase","Isocitrate dehydrogenase (NAD+)","Isocitrate dehydrogenase (NADP+), mitochondrial",
	       #"2-oxoglutarate dehydrogenase complex", "Succinate--CoA ligase (GDP-forming)", "Succinate--CoA ligase (ADP-forming)", "fumarate hydratase",
	       #"malate dehydrogenase 2, NAD (mitochondrial)")


# Electron transport chain complexes I-V
ETC_inds <- c(108:112)
ETC_names <- names[ETC_inds]#c("NADH dehydrogenase", "succinate dehydrogenase", "cytochrome c reductase", "cytochrome c oxidase", "ATP synthase")

# Make a data frame to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
d <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
d2 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
d3 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
d4 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)


# Get a baseline with MitoMammal, because we have a "full" set of mitochondrial functions
mm.base <- optimizeProb(mm, lpdir = "max")
mm.FD <- data.frame(x = getFluxDist(mm.base))
names(mm.FD) <- "Flux (µM/min/gDW)"

# Extract the data
these_data <- data.frame(KO="NONE", EX_O2 = round(mm.FD[EX_o2,],4), MAX_ATP = mm.FD[OBJ_inds[1],], ETC = t(mm.FD[ETC_inds,]), TCA = t(mm.FD[TCA_inds[1],]))

d <- rbind(d, these_data)
d2 <- rbind(d2, these_data)


## Do triple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)
#
#for(i in 1:ncol(KO_inds)){
#  mm.KO <- optimizeProb(mm, react = c(KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(0,0,0), ub = c(0,0,0), lpdir = "max")
#
#  # Fetch the flux distribution
#  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
#  names(mm.KO.FD) <- "Flux (µM/min/gDW)"
#
#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = #t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
#
#  d2 <- rbind(d2, these_data)
#}
#
## Do quadruple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=4)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 4)
#
#for(i in 1:ncol(KO_inds)){
#  mm.KO <- optimizeProb(mm, react = c(KO_inds[1,i], KO_inds[2,i], KO_inds[3,i], KO_inds[4,i]), lb = c(0,0,0,0), ub = c(0,0,0,0), lpdir = "max")
#
#  # Fetch the flux distribution
#  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
#  names(mm.KO.FD) <- "Flux (µM/min/gDW)"
#
#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i], KO_names[4,i], sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = #t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
#
#  d2 <- rbind(d2, these_data)
#}
#
#names(d2) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")
#
#print(d2)
#
## Print to file
#filename <- "triple-quadruple-KO-MAX_ATP-no-restriction.csv"
#write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d2, row.names = F)

## Do triple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)
#
#for(i in 1:ncol(KO_inds)){
#  mm.KO <- optimizeProb(mm, react = c(EX_o2, KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(o20/10,0,0,0), ub = c(-o20/10,0,0,0), lpdir = "max")
#
#  # Fetch the flux distribution
#  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
#  names(mm.KO.FD) <- "Flux (µM/min/gDW)"
#
#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = #t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
#
#  d4 <- rbind(d4, these_data)
#}
#
## Do quadruple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=4)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 4)
#
#for(i in 1:ncol(KO_inds)){
#  mm.KO <- optimizeProb(mm, react = c(EX_o2, KO_inds[1,i], KO_inds[2,i], KO_inds[3,i], KO_inds[4,i]), lb = c(o20/10,0,0,0,0), ub = c(-o20,0,0,0,0), lpdir = "max")
#
#  # Fetch the flux distribution
#  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
#  names(mm.KO.FD) <- "Flux (µM/min/gDW)"
#
#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i], KO_names[4,i], sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = #t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
#
#  d4 <- rbind(d4, these_data)
#}
#
#names(d4) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")
#
#print(d4)
#
## Print to file
#filename <- "triple-quadruple-KO-MAX_ATP-oxygen-restriction.csv"
#write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d4, row.names = F)