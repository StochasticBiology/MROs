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
