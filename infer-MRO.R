#################
source("hypertraps.R")

df.ncbi = read.table("Data/mro-ncbi-tree.phy-cooked.txt-data.txt", header=FALSE)
ancs.ncbi = as.matrix(df.ncbi[seq(1, nrow(df.ncbi), by=2),])
descs.ncbi = as.matrix(df.ncbi[seq(2, nrow(df.ncbi), by=2),])

df.tt = read.table("Data/mro-tree-tt-format.phy-data.txt", header=FALSE)
ancs.tt = as.matrix(df.tt[seq(1, nrow(df.tt), by=2),])
descs.tt = as.matrix(df.tt[seq(2, nrow(df.tt), by=2),])
t1s.tt = as.matrix(read.table("Data/mro-tree-tt-format.phy-datatime.txt", header=FALSE))

df.ttp = read.table("Data/mro-tree-ttplus-format.phy-data.txt", header=FALSE)
ancs.ttp = as.matrix(df.ttp[seq(1, nrow(df.ttp), by=2),])
descs.ttp = as.matrix(df.ttp[seq(2, nrow(df.ttp), by=2),])
t1s.ttp = as.matrix(read.table("Data/mro-ttplus-1.txt", header=FALSE))
t2s.ttp = as.matrix(read.table("Data/mro-ttplus-2.txt", header=FALSE))

features = c("CI", "CII", "CIII", "CIV", "CV", "PDH", "DNA", "TCA", "Fe-S")

### simple demo
my.post.ncbi.1 = HyperTraPS(descs.ncbi, initialstates = ancs.ncbi,
                     losses = 1,
                     length = 5,
                     kernel = 6,
                     seed = 1,
                     samplegap = 10,
                     penalty = 0.,
                     output_transitions = 1,
                     featurenames = features) 
my.post.ncbi.2 = HyperTraPS(descs.ncbi, initialstates = ancs.ncbi,
                          losses = 1,
                          length = 5,
                          kernel = 6,
                          seed = 2,
                          samplegap = 10,
                          penalty = 0.,
                          output_transitions = 1,
                          featurenames = features) 
my.post.ncbi.3 = HyperTraPS(descs.ncbi, initialstates = ancs.ncbi,
                            losses = 1,
                            length = 5,
                            kernel = 6,
                            seed = 3,
                            samplegap = 10,
                            penalty = 0.,
                            output_transitions = 1,
                            featurenames = features) 

sf = 2
png("first-output.png", width=1000*sf, height=800*sf, res=72*sf)
print(ggarrange(plotHypercube.lik.trace(my.post.ncbi.1), 
          plotHypercube.lik.trace(my.post.ncbi.2), 
          plotHypercube.lik.trace(my.post.ncbi.3),
          plotHypercube.sampledgraph2(my.post.ncbi.1, thresh=0.05, use.arc=FALSE, edge.label.size=3, no.times=TRUE) + theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4)),
                    plotHypercube.sampledgraph2(my.post.ncbi.2, thresh=0.05, use.arc=FALSE, edge.label.size=3, no.times=TRUE) + theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4)),
                    plotHypercube.sampledgraph2(my.post.ncbi.3, thresh=0.05, use.arc=FALSE, edge.label.size=3, no.times=TRUE) + theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4)),
        #  plotHypercube.influencegraph(my.post.ncbi.1, cv.thresh = 0.75),
        #  plotHypercube.influencegraph(my.post.ncbi.2, cv.thresh = 0.75),
        #  plotHypercube.influencegraph(my.post.ncbi.3, cv.thresh = 0.75),
        plotHypercube.influences(my.post.ncbi.1, cv.thresh = Inf),
        plotHypercube.influences(my.post.ncbi.2, cv.thresh = Inf),
       plotHypercube.influences(my.post.ncbi.3, cv.thresh = Inf),
          nrow=3, ncol=3))
dev.off()

#### OK up to here
save.image("inferred-7-aug.RData")

# some of these look dodgy -- remove?
if(FALSE) {
plotHypercube.sampledgraph2(my.post.ncbi.1, thresh=0.05, use.arc=FALSE, edge.label.size=3, no.times=TRUE) + 
  theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))

plotHypercube.summary(my.post.ncbi)
plotHypercube.sampledgraph2(my.post.ncbi, thresh=0.05, use.arc=FALSE, edge.label.size=3, no.times=TRUE) + 
  theme(legend.position="none") + expand_limits(x = c(-0.4, 1.4))
plotHypercube.influences(my.post.ncbi, cv.thresh = Inf)
plotHypercube.influencegraph(my.post.ncbi, cv.thresh = 0.4)
plotHypercube.motifseries(my.post.ncbi, c( 0.0001, 0.001, 0.01))
}

my.post.tt.1 = HyperTraPS(descs.tt, initialstates = ancs.tt,
                     starttimes = t1s.tt,
                     losses = 1,
                     length = 5,
                     seed = 1,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = features) 
my.post.tt.2 = HyperTraPS(descs.tt, initialstates = ancs.tt,
                          starttimes = t1s.tt,
                          losses = 1,
                          length = 5,
                          seed = 2,
                          samplegap = 10,
                          output_transitions = 1,
                          featurenames = features) 
my.post.tt.3 = HyperTraPS(descs.tt, initialstates = ancs.tt,
                          starttimes = t1s.tt,
                          losses = 1,
                          length = 5,
                          seed = 3,
                          samplegap = 10,
                          output_transitions = 1,
                          featurenames = features) 

my.post.ttp.1 = HyperTraPS(descs.ttp, initialstates = ancs.ttp, 
                     starttimes = t1s.ttp, endtimes = t2s.ttp, 
                     length = 5,
                     losses = 1,
                     seed = 1,
                     samplegap = 10,
                     output_transitions = 1,
                     featurenames = features)
my.post.ttp.2 = HyperTraPS(descs.ttp, initialstates = ancs.ttp, 
                         starttimes = t1s.ttp, endtimes = t2s.ttp, 
                         length = 5,
                         losses = 1,
                         seed = 2,
                         samplegap = 10,
                         output_transitions = 1,
                         featurenames = features)
my.post.ttp.3 = HyperTraPS(descs.ttp, initialstates = ancs.ttp, 
                         starttimes = t1s.ttp, endtimes = t2s.ttp, 
                         length = 5,
                         losses = 1,
                         seed = 3,
                         samplegap = 10,
                         output_transitions = 1,
                         featurenames = features)

save.image("inferred-7-aug-2.RData")

plotHypercube.summary(my.post.tt)
plotHypercube.summary(my.post.ttp)
