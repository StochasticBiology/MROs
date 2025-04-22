<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
=======
#!/bin/bash
>>>>>>> Stashed changes

# Runs classical HyperTraPS on the full NCBI tree (without timings)
MRO-robert.R mro-ncbi-tree-2025.phy mro-barcodes-2025.csv

# Runs HyperTraPS-CT on the TimeTree tree with and without timings
MRO-robert-CT.R mro-timetree-2025.nwk mro-barcodes-2025.csv

# Runs classical HyperTraPS on the full NCBI tree (without timings) and with 10% of loss patterns associated with observations randomly changed (0s to 1s and 1s to 0s)
MRO-robert-CV.R mro-ncbi-tree-2025.phy mro-barcodes-2025.csv
=======
#!/bin/bash

# Runs classical HyperTraPS on the full NCBI tree (without timings)
MRO-robert.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv &

# Runs HyperTraPS-CT on the TimeTree tree with and without timings
MRO-robert-CT.R mro-timetree-2025.nwk mro-barcodes-2025.csv &

# Produce 10 datasets in which 10% of the loss patterns associated with observations are changed (0s to 1s and 1s to 0s)
randomize-data.R mro-barcodes-2025.csv
=======
#!/bin/bash

# Use: MRO-rcg.R [treefile] [file with binary string] [output folder] [penalty (0/1)] [with_timings (0/1)] 

# Runs classical HyperTraPS on the full NCBI tree (without timings, and without penalty)
MRO-rcg.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv HyperTraPS 0 0 &

# Runs classical HyperTraPS on the full NCBI tree (without timings, with penalty)
MRO-rcg.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv HyperTraPS 1 0 &

# Runs HyperTraPS-CT on the TimeTree tree with and without timings
MRO-rcg.R mro-timetree-2025.nwk mro-barcodes-2025.csv HyperTraPS-CT 0 1 &

# Runs HyperTraPS on the Apicomplexans
MRO-rcg.R mro-tree-apicomplexans.nwk mro-barcodes-2025-apicomplexans.csv Apicomplexans 0 0 &

# Runs HyperTraPS on the non-Apicomplexans
MRO-rcg.R mro-tree-non-apicomplexans.nwk mro-barcodes-2025-non-apicomplexans.csv Non-Apicomplexans 0 0 &

# Makes 10 datasets in which 10% of observations have one entry changed (0s to 1s and 1s to 0s)
randomizeData.R 10 mro-barcodes-2025.csv
>>>>>>> Stashed changes
=======
=======
>>>>>>> Stashed changes
#!/bin/bash

# Use: MRO-rcg.R [treefile] [file with binary string] [output folder] [penalty (0/1)] [with_end_times (0/1)] [with_start_times (0/1)] 

# Runs classical HyperTraPS on the full NCBI tree (without timings, and without penalty) (check: No observations that do not correspond to tips of the tree)
./MRO-rcg.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv HyperTraPS 0 0 &

# Runs classical HyperTraPS on the full NCBI tree (without timings, with penalty) (check)
MRO-rcg.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv HyperTraPS 1 0 &

# Runs HyperTraPS-CT on the TimeTree tree with and without timings (not check?)
./MRO-rcg.R mro-timetree-2025.nwk mro-barcodes-2025.csv HyperTraPS-CT 0 1 &

# Runs HyperTraPS on the Apicomplexans (check)
./MRO-rcg.R mro-ncbi-tree-2025-apicomplexans.nwk mro-barcodes-2025-apicomplexans.csv Apicomplexans 0 0 &

# Runs HyperTraPS on the non-Apicomplexans (check)
./MRO-rcg.R mro-ncbi-tree-2025-non-apicomplexans.nwk mro-barcodes-2025-non-apicomplexans.csv Non-Apicomplexans 0 0 &

# Runs HyperTraPS on the unicellular organisms alone (check)
./MRO-rcg.R mro-ncbi-tree-2025-unicellular.nwk mro-barcodes-2025-unicellular.csv Unicellular 0 0 &

# Runs HyperTraPS on the ciliophora (check)
./MRO-rcg.R mro-ncbi-tree-2025-ciliophora.nwk mro-barcodes-2025-ciliophora.csv Ciliophora 0 0 &
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

# Wait until jobs are finished
wait

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
# Runs classical HyperTraPS on the full NCBI tree (without timings) and with 10% of loss patterns associated with observations randomly changed (0s to 1s and 1s to 0s)
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-1.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-2.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-3.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-4.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-5.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-6.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-7.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-8.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-9.csv &
MRO-robert-CV.R mro-ncbi-tree-2025.nwk mro-barcodes-2025-10.csv &
=======
# Runs HyperTraPS on the ten randomized datasets
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-1.csv HyperTraPS-uncertainty 1 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-2.csv HyperTraPS-uncertainty 2 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-3.csv HyperTraPS-uncertainty 3 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-4.csv HyperTraPS-uncertainty 4 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-5.csv HyperTraPS-uncertainty 5 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-6.csv HyperTraPS-uncertainty 6 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-7.csv HyperTraPS-uncertainty 7 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-8.csv HyperTraPS-uncertainty 8 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-9.csv HyperTraPS-uncertainty 9 &
MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk randomized-dataset-10.csv HyperTraPS-uncertainty 10

>>>>>>> Stashed changes
=======
=======
>>>>>>> Stashed changes
# Runs HyperTraPS on ten randomized datasets

# Use: MRO-rcg-uncertainty.R [treefile] [file with binary string] [output folder]
./MRO-rcg-uncertainty.R mro-ncbi-tree-2025.nwk mro-barcodes-2025.csv HyperTraPS-uncertainty

<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

# Run FBA in parallel with uncertainty analysis above
run-fba.R

# Wait until run-fba.R is finished
wait

# Plot results from FBA
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
plot-fba.R
>>>>>>> Stashed changes
=======
plot-fba.R
>>>>>>> Stashed changes
=======
plot-fba.R
>>>>>>> Stashed changes
=======
plot-fba.R
>>>>>>> Stashed changes
