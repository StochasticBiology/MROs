#!/usr/bin/env Rscript

library(ggplot2)
library(sybilSBML)

fba.mod = readSBMLmod("FBA/MitoMammal/MitoMammal.xml")
etc.inds <- c(108:112)
etc.rs = fba.mod@react_name[etc.inds]

#cat("Reactions are...")
#print(etc.rs)
#cat("END\n\n")

atp.obj = rep(0, length(fba.mod@react_name))
atp.obj[71] = 1
res.df = data.frame(name = NULL, ko.val = NULL, ko.atp.val = NULL, ko.val.hypoxia = NULL, ko.atp.val.hypoxia = NULL)
for(ko in 1:length(fba.mod@react_name)) {
  tmp.mod = fba.mod
  tmp.mod@lowbnd[ko] = 0
  tmp.mod@uppbnd[ko] = 0
  soln = optimizeProb(tmp.mod)
  soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
  tmp.mod@lowbnd[50] = 0
  tmp.mod@uppbnd[50] = 1
  soln.hypoxia = optimizeProb(tmp.mod)
  soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
  res.df = rbind(res.df, data.frame(name=fba.mod@react_name[ko], 
                                    ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                    ko.val.hypoxia=soln.hypoxia@lp_obj,
                                    ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))
}
tmp.mod = fba.mod
soln = optimizeProb(tmp.mod)
soln.atp = optimizeProb(tmp.mod, obj_coef = atp.obj)
tmp.mod@lowbnd[50] = 0
tmp.mod@uppbnd[50] = 1
soln.hypoxia = optimizeProb(tmp.mod)
soln.atp.hypoxia = optimizeProb(tmp.mod, obj_coef = atp.obj)
res.df = rbind(res.df, data.frame(name="full model", 
                                  ko.val=soln@lp_obj, ko.atp.val=soln.atp@lp_obj,
                                  ko.val.hypoxia=soln.hypoxia@lp_obj,
                                  ko.atp.val.hypoxia=soln.atp.hypoxia@lp_obj))


res.df = read.csv("fba-res.csv")
res.df$label = ""
res.df$label[108] = "CI"
res.df$label[109] = "CII"
res.df$label[110] = "CIII"
res.df$label[111] = "CIV"
res.df$label[112] = "CV"

cat("Testing 1...\n")
write.csv(res.df, "fba-res.csv")
ggplot(res.df, aes(x=ko.atp.val, y=ko.atp.val.hypoxia, label=label)) + geom_point() + geom_text()
cat("Testing 2...\n")

ggplot(res.df, aes(x=ko.val, y=ko.val.hypoxia, label=label)) + 
  geom_point() + 
  ggrepel::geom_text_repel()

first.steps = parallelised.runs[[1]]$routes[,1]
first.df = res.df[res.df$label!="",]
print(first.df)
first.df$propn = 0
for(i in 1:nrow(first.df)) {
  ref = which(featurenames == first.df$label[i])-1
  propn = sum(first.steps == ref) / length(first.steps)
  first.df$propn[i] = propn
}
g.normoxic = ggplot(first.df, aes(x=ko.atp.val, y=propn, label=label)) + 
  geom_point() + geom_smooth(method="lm", color="#AAAAFF", fill="#AAAAFF") + 
  ggrepel::geom_text_repel() + theme_light() +
  labs(x = "ATP objective on KO", y = "First loss probability")

g.normoxic
fit.lm = summary(lm(propn~ko.atp.val, data=first.df))

cs = fit.lm$coefficients
tstr = sprintf("Slope %.1e +- %.1e, p=%.2e", cs[2,1], cs[2,2], cs[2,4])

png("mro-first-fba.png", width=300*sf, height=300*sf, res=72*sf)
print(g.normoxic + ggtitle(tstr))
dev.off()