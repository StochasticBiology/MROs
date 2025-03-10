#!/usr/bin/env Rscript

# This is the working directory
require(libSBML)
require(sybilSBML) # Loads sybil
require(rsbml)
require(glpkAPI)
require(Rglpk)
require(combinat)

mc <- readSBMLmod("MitoCore/MitoCore-2017.xml")

# Oxygen: O2t		Boundary conditions - core	o2 transport (diffusion)			M_o2_e --> M_o2_c	lb = 0	ub = 19,8	ATP@MAX_ATP = 19,800
EX_o2 <- 413

OBJ_inds <- c(1:4)
ETC_inds <- c(38:42)
TCA_inds <- c(28:37)

# Make a data frame to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
d <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)


# Get a baseline with MitoMammal across all other cases as well, because we have a "full" set of mitochondrial functions
mc.base <- optimizeProb(mc, lpdir = "max")
mc.FD <- data.frame(x = getFluxDist(mc.base))
names(mc.FD) <- "Flux (µM/min/gDW)"

# Extract the data
these_data <- data.frame(KO="NONE", EX_O2 = round(mc.FD[EX_o2,],4), MAX_ATP = mc.FD[OBJ_inds[1],], ETC = t(mc.FD[ETC_inds,]), TCA = t(mc.FD[TCA_inds[1],]))

d <- rbind(d, these_data)

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c("CI", "CII", "CIII", "CIV", "CV", "PDH")

# Do the single knockouts
for(i in 1:length(KO_inds)){
  mc.KO <- optimizeProb(mc, react = KO_inds[i], lb = 0, ub = 0, lpdir = "max")

  # Fetch the flux distribution
  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
  names(mc.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))

  d <- rbind(d, these_data)
}

# Do double knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m = 2) 
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 2)
for(i in 1:ncol(KO_inds)){
  mc.KO <- optimizeProb(mc, react = c(KO_inds[1,i], KO_inds[2,i]), lb = c(0,0), ub = c(0,0), lpdir = "max")

  # Fetch the flux distribution
  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
  names(mc.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))

  d <- rbind(d, these_data)
}

## Do triple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)

#for(i in 1:ncol(KO_inds)){
#  mc.KO <- optimizeProb(mc, react = c(KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(0,0,0), ub = c(0,0,0), lpdir = "max")
#
#  # Fetch the flux distribution
#  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
#  names(mc.KO.FD) <- "Flux (µM/min/gDW)"
#
#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = #t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))
#
#  d <- rbind(d, these_data)
#}
names(d) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

print(d)

# Print to file
filename <- "MAX_ATP-no-restriction.csv"
write.csv(file = paste("MitoCore/Results/",filename,sep=""), x = d, row.names = F)

######## Repeat with a O2 restricted model
# Make a data frame to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
d2 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)

# Get baseline for anoxic case (10% the oxygen availability required by the basecase (it is negative, hence signs)
o20 <- mc.FD[EX_o2,]
mc.base <- optimizeProb(mc, react = EX_o2, lb = -o20/10, ub = o20/10, lpdir = "max")
mc.FD <- data.frame(x = getFluxDist(mc.base))
names(mc.FD) <- "Flux (µM/min/gDW)"

# Extract the data
these_data <- data.frame(KO="NONE", EX_O2 = round(mc.FD[EX_o2,],4), MAX_ATP = mc.FD[OBJ_inds[1],], ETC = t(mc.FD[ETC_inds,]), TCA = t(mc.FD[TCA_inds[1],]))

d2 <- rbind(d2, these_data)

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c("CI", "CII", "CIII", "CIV", "CV", "PDH")

# Do single knockouts
for(i in 1:length(KO_inds)){
  mc.KO <- optimizeProb(mc, react = c(EX_o2,KO_inds[i]), lb = c(-o20/10,0), ub = c(o20/10,0), lpdir = "max")

  # Fetch the flux distribution
  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
  names(mc.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
}

# Do double knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m = 2) 
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 2)
for(i in 1:ncol(KO_inds)){
  mc.KO <- optimizeProb(mc, react = c(EX_o2,KO_inds[1,i], KO_inds[2,i]), lb = c(-o20/10,0,0), ub = c(o20/10,0,0), lpdir = "max")

  # Fetch the flux distribution
  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
  names(mc.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
}

# Do triple knockouts
#KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
#KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)

#for(i in 1:ncol(KO_inds)){
#  mc.KO <- optimizeProb(mc, react = c(EX_o2, KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(-o20/10,0,0,0), ub = c(o20/10,0,0,0), lpdir = "max")

#  # Fetch the flux distribution
#  mc.KO.FD <- data.frame(x=getFluxDist(mc.KO))
#  names(mc.KO.FD) <- "Flux (µM/min/gDW)"

#  # Extract the data
#  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mc.KO.FD[EX_o2,],4), MAX_ATP = mc.KO.FD[OBJ_inds[1],], ETC = #t(mc.KO.FD[ETC_inds,]), TCA = t(mc.KO.FD[TCA_inds[1],]))

#  d2 <- rbind(d2, these_data)
#}
names(d2) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

print(d2)

# Print to file
filename <- "MAX_ATP-oxygen-restriction.csv"
write.csv(file = paste("MitoCore/Results/",filename,sep=""), x = d2, row.names = F)