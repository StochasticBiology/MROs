# prepare and wrangle data MRO case studies

cd Process
chmod +x *.sh

# prepare MRO-NCBI data
./cook-data.sh ../Data/mro-ncbi-tree.phy ../Data/mro-barcodes.csv 0.1 1
# prepare MRO-TT data (special steps needed)
./mro-timetree-parse.sh
