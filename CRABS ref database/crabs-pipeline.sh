# Create virtual environment for CRABS and dowload all packages required 
conda create -n crabs-env python=3.11.7
conda activate crabs-env
# Install CRABS (v.1.9.0.)
git clone https://github.com/gjeunen/reference_database_creator.git
# Install python packages 
pip install requests
pip install rich
pip install rich-click
pip install matplotlib
pip install numpy

# Install external softwares programs
brew install blast
pip install cutadapt
brew install vsearch
conda install bioconda::clustalw
conda install bioconda::fasttree

# Check the installation
crabs -h

# Before start check if you're in the right paste to save all the files 
# I'm gonna use the MitoFish and NCBI databases

# STEP 1 - DOWNLOAD SEQUENCE DATA FROM ONLINE REPOSITORIES AND TAXONOMY FROM NCBI
# Reference sequences
crabs --download-mitofish --output mitofish.fasta
crabs --download-ncbi --query '("Chondrichthyes"[Organism] OR "Dipnomorpha"[Organism] OR "Actinopterygii"[Organism] OR "Myxini"[Organism] 
OR "Hyperoartia"[Organism] OR "Coelacanthimorpha"[Organism] OR Fish[All Fields]) AND 12S[All Fields]' --output ncbi_fish12s.fasta 
--email eduarda.jesus@usp.br --database nucleotide # need a lot of computer capacity -> I didn't dowload
# Taxonomy from NCBI
crabs --download-taxonomy --exclude 'acc2taxid, taxdump'
crabs --download-taxonomy --output ncbi_taxonomy #ncbi_taxonomy is a subdirectory that I crated -> mkdir ncbi_taxonomy 

# STEP 2 - IMPORT DOWNLOADED DATA INTO CRABS FORMAT
# mitofish
crabs --import --import-format mitofish --input mitofish.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp
--acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_mitofish.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'
# ncbi - If downloaded 
crabs --import --import-format ncbi --input nci.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp
--acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_ncbi.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'

# STEP 3 - MERGE DOWLOADED SEQUENCES
# This step is just when sequence data from multiple online repositories are downloaded. After import the data in crabs format step.
# Not my case
crabs --merge --input 'crabs_mitofish.txt;crabs_ncbi.txt' --uniq --output merged.txt

# STEP 4 - EXTRACT AMPLICONS REGIONS THROUGH IN SILICO PCR ANALYSIS
crabs --in-silico-pcr --input crabs_mitofish.txt --output insilico.txt --forward AAACTCGTGCCAGCCACC --reverse GGGTATCTAATCCCAGTTTG
