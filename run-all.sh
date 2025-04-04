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

# Wait until jobs are finished
wait

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

# Run FBA in parallel with uncertainty analysis above
run-fba.R

# Wait until run-fba.R is finished
wait

# Plot results from FBA
plot-fba.R
>>>>>>> Stashed changes
