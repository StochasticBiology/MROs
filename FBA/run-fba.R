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
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
reacts <- read.csv("MitoMammal/MitoMammal-reactions.txt")

names <- reacts[,2]

# Objective indices
OBJ_inds <- c(71:74)
OBJ_names <- c("OF_ATP_MitoCore","OF_HEME_MitoCore","OF_LIPID_MitoCore","OF_PROTEIN_MitoCore")
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
reacts <- mm@react_name #read.csv("MitoMammal/MitoMammal-reactions.txt")

names <- reacts

# Objective indices
OBJ_inds <- c(71:74)
OBJ_names <- names[OBJ_inds]
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

# Oxygen exchange:
EX_o2 <- 50

# Test that it's really oxygen (check)
# changeObjFunc(mm, react = EX_o2)

# Tricarboxylic acid cycle complexes
TCA_inds <- c(98:107)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
TCA_names <- c("pyruvate dehydrogenase","citrate synthase","Aconitate hydratase","Isocitrate dehydrogenase (NAD+)","Isocitrate dehydrogenase (NADP+), mitochondrial",
	       "2-oxoglutarate dehydrogenase complex", "Succinate--CoA ligase (GDP-forming)", "Succinate--CoA ligase (ADP-forming)", "fumarate hydratase",
	       "malate dehydrogenase 2, NAD (mitochondrial)")


# Electron transport chain complexes I-V
ETC_inds <- c(108:112)
ETC_names <- c("NADH dehydrogenase", "succinate dehydrogenase", "cytochrome c reductase", "cytochrome c oxidase", "ATP synthase")
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
TCA_names <- names[TCA_inds]
=======
TCA_names <- c("PDH","CS","ACONT","ICDHx","ICDHy", "AKGD", "SUCOAS1", "SUCOAS", "FUM", "MDH")#names[TCA_inds]
>>>>>>> Stashed changes
#cat(paste(c(names[TCA_inds],"\n")))

# Electron transport chain complexes I-V
ETC_inds <- c(108:112)
<<<<<<< Updated upstream
ETC_names <- names[ETC_inds]
#cat(paste(c(names[ETC_inds],"\n")))
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

# Make a data frame to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
d <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
d2 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
d3 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)
d4 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, TCA = NULL)


# Get a baseline with MitoMammal across all cases, because we have a "full" set of mitochondrial functions
=======


# Get a baseline with MitoMammal, because we have a "full" set of mitochondrial functions
>>>>>>> Stashed changes
=======


# Get a baseline with MitoMammal, because we have a "full" set of mitochondrial functions
>>>>>>> Stashed changes
=======


# Get a baseline with MitoMammal, because we have a "full" set of mitochondrial functions
>>>>>>> Stashed changes
=======
ETC_names <- c("CI","CII","CIII","CIV","CV")#names[ETC_inds]
#cat(paste(c(names[ETC_inds],"\n")))

# Make a data frame to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
d <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, PDH = NULL)
d2 <- data.frame(KO = NULL, EX_O2 = NULL, MAX_ATP = NULL, ETC = NULL, PDH = NULL)


# Get a baseline with MitoMammal, because we have a "full" set of mitochondrial functions
>>>>>>> Stashed changes
mm.base <- optimizeProb(mm, lpdir = "max")
mm.FD <- data.frame(x = getFluxDist(mm.base))
names(mm.FD) <- "Flux (µM/min/gDW)"

# Extract the data
<<<<<<< Updated upstream
these_data <- data.frame(KO="NONE", EX_O2 = round(mm.FD[EX_o2,],4), MAX_ATP = mm.FD[OBJ_inds[1],], ETC = t(mm.FD[ETC_inds,]), TCA = t(mm.FD[TCA_inds[1],]))

d <- rbind(d, these_data)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
d2 <- rbind(d2, these_data)
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c("CI", "CII", "CIII", "CIV", "CV", "PDH")
=======
these_data <- data.frame(KO="NONE", EX_O2 = round(mm.FD[EX_o2,],4), MAX_ATP = mm.FD[OBJ_inds[1],], ETC = t(mm.FD[ETC_inds,]), PDH = t(mm.FD[TCA_inds[1],]))

d <- rbind(d, these_data)

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c(ETC_names,TCA_names[1])
>>>>>>> Stashed changes

# Do the single knockouts
for(i in 1:length(KO_inds)){
  mm.KO <- optimizeProb(mm, react = KO_inds[i], lb = 0, ub = 0, lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
<<<<<<< Updated upstream
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
=======
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]),PDH = t(mm.KO.FD[TCA_inds[1],]))
>>>>>>> Stashed changes

  d <- rbind(d, these_data)
}

# Do double knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m = 2) 
<<<<<<< Updated upstream
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 2)
=======
KO_names <- combn(c(ETC_names, TCA_names[1]), m = 2)
>>>>>>> Stashed changes
for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(KO_inds[1,i], KO_inds[2,i]), lb = c(0,0), ub = c(0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
<<<<<<< Updated upstream
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))
=======
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), PDH = t(mm.KO.FD[TCA_inds[1],]))
>>>>>>> Stashed changes

  d <- rbind(d, these_data)
}

<<<<<<< Updated upstream
names(d) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
print(d)
=======
#print(d)
>>>>>>> Stashed changes
=======
#print(d)
>>>>>>> Stashed changes
=======
#print(d)
>>>>>>> Stashed changes

# Print to file
filename <- "single-double-KO-MAX_ATP-no-restriction.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d, row.names = F)

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
# Do triple knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)

for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(0,0,0), ub = c(0,0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
}

# Do quadruple knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=4)
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 4)

for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(KO_inds[1,i], KO_inds[2,i], KO_inds[3,i], KO_inds[4,i]), lb = c(0,0,0,0), ub = c(0,0,0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i], KO_names[4,i], sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
}

names(d2) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

print(d2)

# Print to file
filename <- "triple-quadruple-KO-MAX_ATP-no-restriction.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d2, row.names = F)


######## Repeat with a O2 restricted model

# Get baseline for anoxic case (10% the oxygen availability required by the basecase (it is negative, hence signs)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
names(d) <- c("KO","EX_O2","MAX_ATP",ETC_names,TCA_names[1])

#print(d)

# Print to file
filename <- "single-double-KO-MAX_ATP-normoxic.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d, row.names = F)

>>>>>>> Stashed changes

######## Repeat with a O2 restricted model

# Get baseline for hypoxic case (10% the oxygen availability required by the basecase (it is negative, hence signs))
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
o20 <- mm.FD[EX_o2,]
mm.base <- optimizeProb(mm, react = EX_o2, lb = o20/10, ub = -o20/10, lpdir = "max")
mm.FD <- data.frame(x = getFluxDist(mm.base))
names(mm.FD) <- "Flux (µM/min/gDW)"

# Extract the data
<<<<<<< Updated upstream
these_data <- data.frame(KO="NONE", EX_O2 = round(mm.FD[EX_o2,],4), MAX_ATP = mm.FD[OBJ_inds[1],], ETC = t(mm.FD[ETC_inds,]), TCA = t(mm.FD[TCA_inds[1],]))

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
d3 <- rbind(d3, these_data)
d4 <- rbind(d4, these_data)
=======
d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
=======
d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
=======
d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c("CI", "CII", "CIII", "CIV", "CV", "PDH")
=======
these_data <- data.frame(KO="NONE", EX_O2 = round(mm.FD[EX_o2,],4), MAX_ATP = mm.FD[OBJ_inds[1],], ETC = t(mm.FD[ETC_inds,]), PDH = t(mm.FD[TCA_inds[1],]))

d2 <- rbind(d2, these_data)

KO_inds <- c(ETC_inds, TCA_inds[1])
KO_names <- c(ETC_names, TCA_names[1])
>>>>>>> Stashed changes

# Do single knockouts
for(i in 1:length(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(EX_o2,KO_inds[i]), lb = c(o20/10,0), ub = c(-o20/10,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
<<<<<<< Updated upstream
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
  d3 <- rbind(d3, these_data)
=======
  d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
=======
  d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
=======
  d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
=======
  these_data <- data.frame(KO=KO_names[i], EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), PDH = t(mm.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
>>>>>>> Stashed changes
}

# Do double knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m = 2) 
<<<<<<< Updated upstream
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 2)
=======
KO_names <- combn(c(ETC_names, TCA_names[1]), m = 2)
>>>>>>> Stashed changes
for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(EX_o2,KO_inds[1,i], KO_inds[2,i]), lb = c(o20/10,0,0), ub = c(-o20/10,0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
<<<<<<< Updated upstream
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
  d3 <- rbind(d3, these_data)
}

names(d3) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

print(d3)

# Print to file
filename <- "single-double-KO-MAX_ATP-oxygen-restriction.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d3, row.names = F)

# Do triple knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=3)
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 3)

for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(EX_o2, KO_inds[1,i], KO_inds[2,i], KO_inds[3,i]), lb = c(o20/10,0,0,0), ub = c(-o20/10,0,0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

  d4 <- rbind(d4, these_data)
}

# Do quadruple knockouts
KO_inds  <- combn(c(ETC_inds,TCA_inds[1]), m=4)
KO_names <- combn(c("CI", "CII", "CIII", "CIV", "CV", "PDH"), m = 4)

for(i in 1:ncol(KO_inds)){
  mm.KO <- optimizeProb(mm, react = c(EX_o2, KO_inds[1,i], KO_inds[2,i], KO_inds[3,i], KO_inds[4,i]), lb = c(o20/10,0,0,0,0), ub = c(-o20,0,0,0,0), lpdir = "max")

  # Fetch the flux distribution
  mm.KO.FD <- data.frame(x=getFluxDist(mm.KO))
  names(mm.KO.FD) <- "Flux (µM/min/gDW)"

  # Extract the data
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i], KO_names[3,i], KO_names[4,i], sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), TCA = t(mm.KO.FD[TCA_inds[1],]))

  d4 <- rbind(d4, these_data)
}

names(d4) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")

print(d4)

# Print to file
filename <- "triple-quadruple-KO-MAX_ATP-oxygen-restriction.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d4, row.names = F)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
  d2 <- rbind(d2, these_data)
}

names(d2) <- c("KO","EX_O2","MAX_ATP","CI","CII","CIII","CIV","CV","PDH")
=======
  these_data <- data.frame(KO=paste(KO_names[1,i], KO_names[2,i],sep=","), EX_O2 = round(mm.KO.FD[EX_o2,],4), MAX_ATP = mm.KO.FD[OBJ_inds[1],], ETC = t(mm.KO.FD[ETC_inds,]), PDH = t(mm.KO.FD[TCA_inds[1],]))

  d2 <- rbind(d2, these_data)
}

names(d2) <- c("KO","EX_O2","MAX_ATP",ETC_names,TCA_names[1])
>>>>>>> Stashed changes

#print(d2)

# Print to file
<<<<<<< Updated upstream
filename <- "single-double-KO-MAX_ATP-oxygen-restriction.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d2, row.names = F)

# Make a heatmap showing the complexes, with (NONE, KOs) on the y-axis and values of (NONE, KOs) after KO for MAX_ATP in normoxic and anoxic conditions
d.norm <- d[1:6,c(1,3:8)]
d.anox <- d2[1:6,c(1,3:8)]

for(i in 2:6){
d.anox[,i]<-d.anox[,i]/d.norm[1,i]
d.norm[,i]<-d.norm[,i]/d.norm[1,i]
}

print(d.norm)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
print(d.anox)
>>>>>>> Stashed changes
=======
print(d.anox)
>>>>>>> Stashed changes
=======
print(d.anox)
>>>>>>> Stashed changes
=======
filename <- "single-double-KO-MAX_ATP-anoxic.csv"
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = d2, row.names = F)
>>>>>>> Stashed changes
