# MROs

Inference and investigation of reductive evolution of mitochondrion-related organelles (MROs)

Broadly, the pipeline:
* Acquires and curates data on MRO features and their phylogenetic embedding
* Runs hypercubic inference on this data to explore evolutionary pathways
* Looks at metabolic predictors of these pathways using flux balance analysis

At the moment, `MRO-all-R-mac.R` does part of all of these. Given source data (in `Data/`) on an estimated (no branch lengths) phylogeny and presence/absence of MRO features, it uses functionality from HyperTraPS-CT (https://github.com/StochasticBiology/hypertraps-ct) to reconstruct a set of evolutionary transitions. It then runs several parallelised instances of HyperTraPS to look at evolutionary dynamics. Various visualisations are produced.

We hypothesise a link between metabolic influence of a feature and its propensity to be lost. Due to architectural awkwardness in getting linear programming solvers working, `MRO-all-R-mac.R` doesn't do FBA itself, instead just retrieving previous results. This needs improving and at least some self-contained FBA implementation should be included in the pipeline.

Outside `MRO-all-R-mac.R` there are various other files in Bash, Python, etc that form part of the older pipeline for this analysis. There are several choices in the pipeline -- branch lengths vs not, which source of phylogenetic information to use, etc. Ideally all this will get folded in to the core R file and these older files will become deprecated.
