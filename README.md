# MROs

Inference and investigation of reductive evolution of mitochondrion-related organelles (MROs).

<a href="https://doi.org/10.5281/zenodo.16798558"><img src="https://zenodo.org/badge/826822041.svg" alt="DOI"></a>

![image](https://github.com/user-attachments/assets/37c59f21-8094-4c48-b8a9-2c60a3b6a3a5)

Dependencies
----
This code requires R, sybil, sybilSBML, and HyperTraPS-CT and its dependencies (see https://github.com/StochasticBiology/hypertraps-ct). The code is, by default, run in parallel; with 8 cores it should take perhaps 20-30min for the whole pipeline. 

Outline
----
Broadly, the pipeline:
* Curates data on MRO features and their phylogenetic embedding
* Runs hypercubic inference on this data to explore evolutionary pathways
* Looks at metabolic predictors of these pathways using flux balance analysis

`prepare-pub.sh` is the wrapper script for the pipeline. Given source data (in `Data/`) on presence/absence of MRO features and phylogenies with and without branch lengths, it uses functionality from HyperTraPS-CT (https://github.com/StochasticBiology/hypertraps-ct) to reconstruct a set of evolutionary transitions. It examines dynamics in different subsets of the data, reflecting different clades, and also creates several resampled datasets to act as controls, and runs inference on these. Various visualisations are produced.

We hypothesise a link between metabolic influence of a feature and its propensity to be lost. The associated flux balance analysis is performed by the script `FBA/run-fba.R` (also called from `prepare-pub.sh`). 

Outside `prepare-pub.R`, `summarizeHyperTraPS-inference.R` outputs visualizations with which to determined if chains mix well, if output converges, among other things. 

Data
----
Data with references and commentary can be found in `MRO-metaanalysis-2025.xlsx`. This is intepreted according to the various binarisation protocols described therein and (optionally) subsetted by clade, and the binarised curated results are stored in `Data/` for the analysis pipeline.
