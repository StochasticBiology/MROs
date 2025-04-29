# MROs

Inference and investigation of reductive evolution of mitochondrion-related organelles (MROs)

Broadly, the pipeline:
* Acquires and curates data on MRO features and their phylogenetic embedding
* Runs hypercubic inference on this data to explore evolutionary pathways
* Looks at metabolic predictors of these pathways using flux balance analysis

At the moment, `prepare-pub.sh` does part of all of these. Given source data (in `Data/`) on an estimated (no branch lengths) phylogeny and presence/absence of MRO features, it uses functionality from HyperTraPS-CT (https://github.com/StochasticBiology/hypertraps-ct) to reconstruct a set of evolutionary transitions. It then runs several parallelised instances of HyperTraPS to look at evolutionary dynamics. Various visualisations are produced.

We hypothesise a link between metabolic influence of a feature and its propensity to be lost. `prepare-pub.sh` takes a first look at this with a simple model-solver combo in the script `summarizeHyperTraPS-inference.R`. A broader flux balance analysis is performed by running the script `FBA/run-fba.R` (will fold into overall pipeline). 

Outside `prepare-pub.R`, `summarizeHyperTraPS-inference.R` outputs visualizations with which to determined if chains mix well, if output converges, among other things. (Deprecated: There are various other files in Bash, Python, etc that form part of the older pipeline for this analysis. There are several choices in the pipeline -- branch lengths vs not, which source of phylogenetic information to use, etc. Ideally all this will get folded in to the core R file and these older files will become deprecated.

`MRO-all-R-mac.R` reads the FBA work from a premade file, as I can't get the FBA solver working on a Mac.)
