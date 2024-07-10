#################
source("hypertraps.R")

df.ncbi = read.table("Data/mro-ncbi-tree.phy-cooked.txt-data.txt", header=FALSE)
ancs.ncbi = as.matrix(df[seq(1, nrow(df.ncbi), by=2),])
descs.ncbi = as.matrix(df[seq(2, nrow(df.ncbi), by=2),])

df.tt = read.table("Data/mro-tree-tt-format.phy-data.txt", header=FALSE)
ancs.tt = as.matrix(df[seq(1, nrow(df.tt), by=2),])
descs.tt = as.matrix(df[seq(2, nrow(df.tt), by=2),])
t1s.tt = as.matrix(read.table("Data/mro-tree-tt-format.phy-datatime.txt", header=FALSE))

df.ttp = read.table("Data/mro-tree-ttplus-format.phy-data.txt", header=FALSE)
ancs.ttp = as.matrix(df[seq(1, nrow(df.ttp), by=2),])
descs.ttp = as.matrix(df[seq(2, nrow(df.ttp), by=2),])
t1s.ttp = as.matrix(read.table("Data/mro-ttplus-1.txt", header=FALSE))
t2s.ttp = as.matrix(read.table("Data/mro-ttplus-2.txt", header=FALSE))

features = c("CI", "CII", "CIII", "CIV", "CV", "PDH", "DNA", "TCA", "Fe-S")

### simple demo
my.post.ncbi = HyperTraPS(descs.ncbi, initialstates = ancs.ncbi,
                     losses = 1,
                     length = 4,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = features) 

plotHypercube.summary(my.post)
plotHypercube.sampledgraph2(my.post, thresh=0.1, use.arc=FALSE, edge.label.size=3) + 
  theme(legend.position="none") + expand_limits(x = c(-0.1, 1.1))
plotHypercube.influences(my.post, cv.thresh = Inf)
plotHypercube.influencegraph(my.post, cv.thresh = 0.5)
plotHypercube.motifseries(my.post, c( 0.0001, 0.001, 0.01))

my.post.tt = HyperTraPS(descs.tt, initialstates = ancs.tt,
                     starttimes = t1s.tt,
                     losses = 1,
                     length = 4,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = features) 

my.post.ttp = HyperTraPS(descs.ttp, initialstates = ancs.ttp, 
                     starttimes = t1s.ttp, endtimes = t2s.ttp, 
                     length = 4,
                     losses = 1,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = c("A", "B", "C", "D", "E")); 

