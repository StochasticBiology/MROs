#!/bin/bash

# make everything executable
chmod +x *.R
chmod +x FBA/*.R

# get the HyperTraPS-CT source code, and put the essential code into this directory
git clone https://github.com/StochasticBiology/hypertraps-ct
cp hypertraps-ct/hypertraps* .

# Use: runHyperTraPS-inference.R [treefile] [file with binary string] [output folder] [penalty (0/1)] [with_timings (0/1)]

#### Run all types of inferences with the interpretation that 1b --> 1

# Runs classical HyperTraPS on the full NCBI tree (without timings, and without penalty)
./runHyperTraPS-inference.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-2.csv HyperTraPS-2 0 0 &

# Runs HyperTraPS-CT on the TimeTree tree with and without timings with two different intervals (all three stored in the same RData file)
./runHyperTraPS-inference.R mro-timetree-2025.nwk mro-barcodes-CT-2025-2.csv HyperTraPS-CT-2 0 1 &

# Runs HyperTraPS on the Apicomplexans
./runHyperTraPS-inference.R mro-ncbi-tree-2025-apicomplexans.nwk mro-barcodes-2025-apicomplexans-2.csv Apicomplexans-2 0 0 &

# Runs HyperTraPS on the ciliophora
./runHyperTraPS-inference.R mro-ncbi-tree-2025-ciliophora.nwk mro-barcodes-2025-ciliophora-2.csv Ciliophora-2 0 0 &

# Make ten datasets with one observation changed
# Use: randomizeData.R [integer number of datasets] [file with binary strings] [label (1: 1b --> 0; 2: 1b --> 1)]
./randomizeData.R 10 mro-barcodes-2025-2.csv 2

# Wait until jobs are finished
wait

# Runs HyperTraPS on the ten randomized datasets
# Use: MRO-rcg-uncertainty.R [treefile] [file with binary string] [output folder] [label]
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-1.csv HyperTraPS-uncertainty-2 1 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-2.csv HyperTraPS-uncertainty-2 2 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-3.csv HyperTraPS-uncertainty-2 3 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-4.csv HyperTraPS-uncertainty-2 4 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-5.csv HyperTraPS-uncertainty-2 5 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-6.csv HyperTraPS-uncertainty-2 6 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-7.csv HyperTraPS-uncertainty-2 7 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-8.csv HyperTraPS-uncertainty-2 8 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-9.csv HyperTraPS-uncertainty-2 9 &
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-10.csv HyperTraPS-uncertainty-2 10

wait

cd FBA
./run-fba.R MAX_ATP AOX-NDH2

wait

cd ..

# Plot results from HyperTraPS uncertainty analysis and produce other figures (still needs the fba figures from FBA/MitoMammal/Results/MAX_ATP/single-double-AOX-NDH2.csv, but everything else should work now, as well as the missing figures fig-4, SI-fig-3, and SI-fig-5)
./prepare-pub.R 2 .7
