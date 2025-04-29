#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

# test if there's exactly two arguments; if not, return an error
if (length(args)!=2) {
  stop("Need two arguments: an output folder, and an output label.\n Usage: run-fba.R [output folder] [output label]\n", call.=FALSE)
}

# This is the working directory
require(libSBML)
require(sybilSBML) # Loads sybil
require(rsbml)
require(glpkAPI)
require(Rglpk)
require(combinat)

obj.name = as.character(args[1])
outfolder = as.character(args[1])
outfolder = paste(outfolder,"/",sep="")
outlabel = as.character(args[2])

# If output label is empty, AOX, NDH2 won't be activated
if(outlabel != ""){
  outlabel = paste("-", outlabel, sep = "")
  # If output label specifies NDH2 and/or AOX are to be active, set their upper bounds to 1000
  if(grepl(x = outlabel, "NDH2", fixed = T)){
    ndh2.lowbnd <- 0
    ndh2.uppbnd <- 1000
  }else{
    ndh2.lowbnd <- 0
    ndh2.uppbnd <- 0
  }
  if(grepl(x = outlabel, "AOX", fixed = T)){
    aox.lowbnd <- 0
    aox.uppbnd <- 1000
  }else{
    aox.lowbnd <- 0
    aox.uppbnd <- 0
  }
}else{
  ndh2.lowbnd <- 0
  ndh2.uppbnd <- 0
  aox.lowbnd <- 0
  aox.uppbnd <- 0
}

# Load up the MitoMammal GEM
fba.mod <- readSBMLmod("MitoMammal/MitoMammal.xml")

# Add AOX, NDH2, setting flux lb=ub=0
fba.mod <- addReact(fba.mod, id = "AOX", met = c("q10h2[m]","o2[m]","q10[m]","h2o[m]"), Scoef = c(-2,-1,2,2), reversible = F, lb = 0, ub = 0)
fba.mod <- addReact(fba.mod, id = "NDH2", met = c("h[m]","nadh[m]","q10[m]","nad[m]","q10h2[m]"),
                    Scoef = c(-1,-1,-.999,1,.999), reversible = F, lb = 0, ub = 0)

obj = rep(0,length(fba.mod@react_name))
if(obj.name == "MAX_ATP"){
  obj.ind = 71
}else if(obj.name == "MAX_LIPID"){
  obj.ind = 73
}else if(obj.name == "MAX_PROTEIN"){
  obj.ind = 74
}else{
  stop("Need input label: MAX_ATP, MAX_LIPID, or MAX_PROTEIN", call. = FALSE)
}
fba.mod = changeObjFunc(fba.mod, react = obj.ind)
obj[obj.ind] = 1

# Oxygen exchange:
EX_o2 <- 50

# AOX and NDH2
AOX_ind =  561
NDH2_ind =  562

# Test that it's really oxygen (check)
# changeObjFunc(fba.mod, react = EX_o2)

# Tricarboxylic acid cycle complexes
TCA_inds <- c(98:107)
TCA_names <- c("PDH","CS","ACONT","ICDHx","ICDHy", "AKGD", "SUCOAS1", "SUCOAS", "FUM", "MDH")

# Electron transport chain complexes I-V
ETC_inds <- c(108:112)
ETC_names <- c("CI","CII","CIII","CIV","CV")

# Make data frames to keep track of what's KOed, and the fluxes at each complex of (primary) interest (thus far, at least)
df <- data.frame(KO = NULL, OBJ = NULL,
			  max_obj_normoxic =  NULL, max_obj_normoxic2 = NULL,
			  max_obj_hypoxic = NULL, max_obj_hypoxic2 = NULL)

# Keep track of fluxes through introduced reactions
df.alts <- data.frame(KO = NULL, OBJ = NULL,
                      AOX.normoxic = NULL, AOX.hypoxic = NULL,
                      NDH2.normoxic = NULL, NDH2.hypoxic = NULL,
                      ETC.normoxic = NULL, ETC.hypoxic = NULL)

# Get a baseline with the full set of mitochondrial functions and save the oxygen exchange for optimal conditions
tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
f.o2.u <- -getFluxDist(soln)[EX_o2]/10
f.o2.l <- -f.o2.u/10
#print(f.o2.u)
# Use 10% of oxygen exchange in optimal conditions for anoxic conditions
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
df = rbind(df, data.frame(KO="None", OBJ = obj.name,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO = "None", OBJ = obj.name,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Setup for single KOs of the ETCs
KOs <- data.frame(ind = c(ETC_inds),
                  name = c(ETC_names))

# Do the single knockouts

for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KO
  KOind = KOs$ind[ko]
  KOname = KOs$name[ko]
  tmp.mod@lowbnd[KOind] = 0
  tmp.mod@uppbnd[KOind] = 0
  
  # If CI is KOed, add NDH2
  if(KOname == "CI"){
   #print(KOname)
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII (the fourth KO) or CIV (the fifth), add AOX
  if(KOname == "CIII" | KOname == "CIV"){
   #print(KOname)
   tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
   tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
  df = rbind(df, data.frame(KO=KOname, OBJ = obj.name,
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO = KOname, OBJ = obj.name,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add single TCA KOs to single KOs
for(ko in 1:length(TCA_inds)){
  tmp.mod = fba.mod
  tmp.mod@lowbnd[TCA_inds[ko]] = 0
  tmp.mod@uppbnd[TCA_inds[ko]] = 0
  KOname = TCA_names[ko]
  soln = optimizeProb(tmp.mod)
  soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
  df = rbind(df, data.frame(KO=paste("Single TCA", KOname, sep=":"), OBJ = obj.name, 
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste("Single TCA", KOname, sep=":"), OBJ = obj.name,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add full TCA to single KOs (Note: TCA includes PDH in MitoMammal's annotation)
tmp.mod = fba.mod
tmp.mod@lowbnd[TCA_inds] = 0
tmp.mod@uppbnd[TCA_inds] = 0
soln = optimizeProb(tmp.mod)
soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
df = rbind(df, data.frame(KO="Full TCA", OBJ = obj.name,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="Full TCA", OBJ = obj.name,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Do double knockouts
KOs <- data.frame(ind1 = c(rep(ETC_inds[1],5), rep(ETC_inds[2],4), rep(ETC_inds[3],3), rep(ETC_inds[4],2), ETC_inds[5]),
                  ind2 = c(c(TCA_inds[1], ETC_inds[2:5]),
                           c(TCA_inds[1], ETC_inds[3:5]),
                           c(TCA_inds[1], ETC_inds[4:5]),
                           c(TCA_inds[1], ETC_inds[5]),
                           TCA_inds[1]),
                  name1 = c(rep(ETC_names[1],5), rep(ETC_names[2],4), rep(ETC_names[3],3), rep(ETC_names[4],2), ETC_names[5]),
                  name2 = c(c(TCA_names[1], ETC_names[2:5]),
                            c(TCA_names[1], ETC_names[3:5]),
                            c(TCA_names[1], ETC_names[4:5]),
                            c(TCA_names[1], ETC_names[5]),
                            TCA_names[1]))
for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KOs
  KOind1 = KOs$ind1[ko]
  KOind2 = KOs$ind2[ko]
  KOname1 = KOs$name1[ko]
  KOname2 = KOs$name2[ko]
  tmp.mod@lowbnd[c(KOind1,KOind2)] = 0
  tmp.mod@uppbnd[c(KOind1,KOind2)] = 0
  
  # If CI is KOed, add NDH2
  if(KOname1 == "CI" | KOname2 == "CI"){
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII or CIV is KOed, add AOX
  if(KOname1 == "CIII" | KOname1 == "CIV" | KOname2 == "CIII" | KOname2 == "CIV"){
    tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
    tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  soln = optimizeProb(tmp.mod)
  soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
  df = rbind(df, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name,
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# KO each ETC with TCA
KOs <- data.frame(ind = ETC_inds,
                  name1 = ETC_names,
                  name2 = rep("TCA",5))
for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KOs
  KOind1 = KOs$ind[ko]
  KOind2 = TCA_inds
  KOname1 = KOs$name1[ko]
  KOname2 = "TCA"
  tmp.mod@lowbnd[c(KOind1,KOind2)] = 0
  tmp.mod@uppbnd[c(KOind1,KOind2)] = 0
  
  # If CI is KOed, add NDH2
  if(KOname1 == "CI"){
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII or CIV is KOed, add AOX
  if(KOname1 == "CIII" | KOname1 == "CIV"){
    tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
    tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
  df = rbind(df, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name, 
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add {CI,CII,CV}, {CIII,CIV,CV}, and entire ETC
tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds[c(1,3:4)]] = 0
tmp.mod@uppbnd[ETC_inds[c(1,3:4)]] = 0
soln = optimizeProb(tmp.mod)
soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
df = rbind(df, data.frame(KO="CI,CII,CV", OBJ = obj.name,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="CI,CII,CV", OBJ = obj.name,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds[c(1:2,5)]] = 0
tmp.mod@uppbnd[ETC_inds[c(1:2,5)]] = 0
soln = optimizeProb(tmp.mod)
soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
df = rbind(df, data.frame(KO="CI,CIII,CIV", OBJ = obj.name,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="CI,CIII,CIV", OBJ = obj.name,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds] = 0
tmp.mod@uppbnd[ETC_inds] = 0
soln = optimizeProb(tmp.mod)
soln.obj = optimizeProb(tmp.mod, obj_coef = obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.obj.hypoxia = optimizeProb(tmp.mod, obj_coef = obj)
df = rbind(df, data.frame(KO="Full ETC", OBJ = obj.name,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.obj@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.obj.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="Full ETC", OBJ = obj.name,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Set up for MAX PMF AND OBJ
obj.name.pmf = "MAX_PMF"
pmf.obj = rep(0,length(fba.mod@react_name))
pmf.obj[c(108,110,111)] = 4
pmf.obj[112] = -2.97

fba.mod <- changeObjFunc(fba.mod, react = c(obj.ind,108,110,111,112), obj_coef = c(1,4,4,4,-2.97))

# Get a baseline with the full set of mitochondrial functions and save the oxygen exchange for optimal conditions
tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
f.o2.u <- -getFluxDist(soln)[EX_o2]/10
f.o2.l = -f.o2.u/10
# Use 10% of oxygen exchange in optimal conditions for anoxic conditions
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
df = rbind(df, data.frame(KO="None", OBJ = obj.name.pmf,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO = KOname, OBJ = obj.name.pmf,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Setup for single KOs of the ETCs
KOs <- data.frame(ind = c(ETC_inds),
                  name = c(ETC_names))

# Do the single knockouts

for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KO
  KOind = KOs$ind[ko]
  KOname = KOs$name[ko]
  tmp.mod@lowbnd[KOind] = 0
  tmp.mod@uppbnd[KOind] = 0
  
  # If CI is KOed, add NDH2
  if(KOname == "CI"){
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII (the fourth KO) or CIV (the fifth), add AOX
  if(KOname == "CIII" | KOname == "CIV"){
    tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
    tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  df = rbind(df, data.frame(KO=KOname, OBJ = obj.name.pmf,
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO = KOname, OBJ = obj.name.pmf,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add single TCA KOs to single KOs
for(ko in 1:length(TCA_inds)){
  tmp.mod = fba.mod
  tmp.mod@lowbnd[TCA_inds[ko]] = 0
  tmp.mod@uppbnd[TCA_inds[ko]] = 0
  KOname = TCA_names[ko]
  
  soln = optimizeProb(tmp.mod)
  soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  df = rbind(df, data.frame(KO=paste("Single TCA", KOname, sep=":"), OBJ = obj.name.pmf,
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste("Single TCA", KOname, sep=":"), OBJ = obj.name.pmf,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add full TCA to single KOs (Note: TCA includes PDH in MitoMammal's annotation)
tmp.mod = fba.mod
tmp.mod@lowbnd[TCA_inds] = 0
tmp.mod@uppbnd[TCA_inds] = 0
soln = optimizeProb(tmp.mod)
soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
df = rbind(df, data.frame(KO="Full TCA", OBJ = obj.name.pmf,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="Full TCA", OBJ = obj.name.pmf,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Do double knockouts
KOs <- data.frame(ind1 = c(rep(ETC_inds[1],5), rep(ETC_inds[2],4), rep(ETC_inds[3],3), rep(ETC_inds[4],2), ETC_inds[5]),
                  ind2 = c(c(TCA_inds[1], ETC_inds[2:5]),
                           c(TCA_inds[1], ETC_inds[3:5]),
                           c(TCA_inds[1], ETC_inds[4:5]),
                           c(TCA_inds[1], ETC_inds[5]),
                           TCA_inds[1]),
                  name1 = c(rep(ETC_names[1],5), rep(ETC_names[2],4), rep(ETC_names[3],3), rep(ETC_names[4],2), ETC_names[5]),
                  name2 = c(c(TCA_names[1], ETC_names[2:5]),
                            c(TCA_names[1], ETC_names[3:5]),
                            c(TCA_names[1], ETC_names[4:5]),
                            c(TCA_names[1], ETC_names[5]),
                            TCA_names[1]))
for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KOs
  KOind1 = KOs$ind1[ko]
  KOind2 = KOs$ind2[ko]
  KOname1 = KOs$name1[ko]
  KOname2 = KOs$name2[ko]
  tmp.mod@lowbnd[c(KOind1,KOind2)] = 0
  tmp.mod@uppbnd[c(KOind1,KOind2)] = 0
  
  # If CI is KOed, add NDH2
  if(KOname1 == "CI" | KOname2 == "CI"){
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII or CIV is KOed, add AOX
  if(KOname1 == "CIII" | KOname1 == "CIV" | KOname2 == "CIII" | KOname2 == "CIV"){
    tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
    tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  df = rbind(df, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name.pmf, 
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name.pmf, 
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# KO each ETC with TCA
KOs <- data.frame(ind = ETC_inds,
                  name1 = ETC_names,
                  name2 = rep("TCA",5))
for(ko in 1:nrow(KOs)){
  tmp.mod = fba.mod
  
  # Do the KOs
  KOind1 = KOs$ind[ko]
  KOind2 = TCA_inds
  KOname1 = KOs$name1[ko]
  KOname2 = "TCA"
  tmp.mod@lowbnd[c(KOind1,KOind2)] = 0
  tmp.mod@uppbnd[c(KOind1,KOind2)] = 0
  
  # If CI is KOed, add NDH2
  if(KOname1 == "CI"){
   tmp.mod@lowbnd[NDH2_ind] = ndh2.lowbnd
   tmp.mod@uppbnd[NDH2_ind] = ndh2.uppbnd
  }
  
  # If CIII or CIV is KOed, add AOX
  if(KOname1 == "CIII" | KOname1 == "CIV"){
    tmp.mod@lowbnd[AOX_ind] = aox.lowbnd
    tmp.mod@uppbnd[AOX_ind] = aox.uppbnd
  }
  
  soln = optimizeProb(tmp.mod)
  soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  tmp.mod@lowbnd[EX_o2] = f.o2.l
  tmp.mod@uppbnd[EX_o2] = f.o2.u
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
  df = rbind(df, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name.pmf,
			    max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			    max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))
  
  fd.normoxic <- getFluxDist(soln)
  fd.hypoxic  <- getFluxDist(soln.hypoxia)
  df.alts <- rbind(df.alts, data.frame(KO=paste(KOname1, KOname2, sep = ","), OBJ = obj.name.pmf,
                                       AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                       NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                       ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))
}

# Add {CI,CII,CV}, {CIII,CIV,CV}, and entire ETC
tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds[c(1,3:4)]] = 0
tmp.mod@uppbnd[ETC_inds[c(1,3:4)]] = 0
soln = optimizeProb(tmp.mod)
soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
df = rbind(df, data.frame(KO="CI,CII,CV", OBJ = obj.name.pmf,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="CI,CII,CV", OBJ = obj.name.pmf,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds[c(1:2,5)]] = 0
tmp.mod@uppbnd[ETC_inds[c(1:2,5)]] = 0
soln = optimizeProb(tmp.mod)
soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.opmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
df = rbind(df, data.frame(KO="CI,CIII,CIV", OBJ = obj.name.pmf,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="CI,CIII,CIV", OBJ = obj.name.pmf,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

tmp.mod = fba.mod
tmp.mod@lowbnd[ETC_inds] = 0
tmp.mod@uppbnd[ETC_inds] = 0
soln = optimizeProb(tmp.mod)
soln.pmf = optimizeProb(tmp.mod, obj_coef = pmf.obj)
tmp.mod@lowbnd[EX_o2] = f.o2.l
tmp.mod@uppbnd[EX_o2] = f.o2.u
soln.hypoxia = optimizeProb(tmp.mod)
soln.pmf.hypoxia = optimizeProb(tmp.mod, obj_coef = pmf.obj)
df = rbind(df, data.frame(KO="Full ETC", OBJ = obj.name.pmf,
			  max_obj_normoxic = soln@lp_obj, max_obj_normoxic2 = soln.pmf@lp_obj,
  			  max_obj_hypoxic = soln.hypoxia@lp_obj, max_obj_hypoxic2 = soln.pmf.hypoxia@lp_obj))

fd.normoxic <- getFluxDist(soln)
fd.hypoxic  <- getFluxDist(soln.hypoxia)
df.alts <- rbind(df.alts, data.frame(KO="Full ETC", OBJ = obj.name.pmf,
                                     AOX.normoxic = fd.normoxic[AOX_ind], AOX.hypoxic = fd.hypoxic[AOX_ind],
                                     NDH2.normoxic = fd.normoxic[NDH2_ind], NDH2.hypoxic = fd.hypoxic[NDH2_ind],
                                     ETC.normoxic = t(fd.normoxic[ETC_inds]), ETC.hypoxic = t(fd.hypoxic[ETC_inds])))

# Print df to file
if(!dir.exists("MitoMammal/Results"))dir.create("MitoMammal/Results")
if(!dir.exists(paste("MitoMammal/Results/",outfolder, sep = "")))dir.create(paste("MitoMammal/Results/",outfolder, sep = ""))
filename <- paste(outfolder, "single-double", outlabel, ".csv", sep = "")
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = df, row.names = F)

names(df.alts) <- c("KO","OBJ",
                    paste("AOX",c("norm","hyp"),sep="."),
                    paste("NDH2",c("norm","hyp"),sep="."),
                    paste(ETC_names,"norm",sep="."),
                    paste(ETC_names,"hyp",sep="."))

df[,c(3:ncol(df))] <- round(df[,c(3:ncol(df))], 2)
df.alts[,c(3:ncol(df.alts))] <- round(df.alts[,c(3:ncol(df.alts))], 2) 

# Print df.alts to file
filename <- paste(outfolder, "single-double", outlabel, "-alts.csv", sep = "")
write.csv(file = paste("MitoMammal/Results/",filename,sep=""), x = df.alts, row.names = F)
