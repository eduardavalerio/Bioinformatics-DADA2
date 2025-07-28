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
crabs -h # In my case shown "command not found" so I needed to add the crabs script in my PATH 

# HOW TO ADD CRABS IN YOUR PATH
ls -l crabs.py  # Look for crabs.py
chmod +x crabs  # Make it executable
./crabs -h # Run directly, if it works you can add in your path permanently
echo $(pwd)
export PATH="$PATH:/Users/eduardavalerio/reference_database_creator" # Replace /Users/eduarda...
source ~/.bashrc # Reload the shell
crabs -h # Should now work from anywhere


################################################################################################################################################
# Before start check if you're in the right paste to save all the files 
# I'm gonna use the MitoFish and NCBI databases

# STEP 1 - DOWNLOAD SEQUENCE DATA FROM ONLINE REPOSITORIES AND TAXONOMY FROM NCBI
# Reference sequences
crabs --download-mitofish --output mitofish.fasta
crabs --download-ncbi --query '("Chondrichthyes"[Organism] OR "Dipnomorpha"[Organism] OR "Actinopterygii"[Organism] OR "Myxini"[Organism] OR "Hyperoartia"[Organism] OR "Coelacanthimorpha"[Organism] OR Fish[All Fields]) AND 12S[All Fields]' --output ncbi_fish12s.fasta --email eduarda.jesus@usp.br --database nucleotide # need a lot of computer capacity -> I didn't dowload
# Taxonomy from NCBI
crabs --download-taxonomy --exclude 'acc2taxid, taxdump'
crabs --download-taxonomy --output ncbi_taxonomy #ncbi_taxonomy is a subdirectory that I crated -> mkdir ncbi_taxonomy 

# STEP 2 - IMPORT DOWNLOADED DATA INTO CRABS FORMAT
# mitofish
crabs --import --import-format mitofish --input mitofish.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp --acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_mitofish.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'
# ncbi - If downloaded 
crabs --import --import-format ncbi --input nci.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp --acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_ncbi.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'

# STEP 3 - MERGE DOWLOADED SEQUENCES
# This step is just when sequence data from multiple online repositories are downloaded. After import the data in crabs format step.
# Not my case
crabs --merge --input 'crabs_mitofish.txt;crabs_ncbi.txt' --uniq --output merged.txt

# STEP 4 - EXTRACT AMPLICONS REGIONS THROUGH IN SILICO PCR ANALYSIS
crabs --in-silico-pcr --input crabs_mitofish.txt --output insilico.txt --forward AAACTCGTGCCAGCCACC --reverse GGGTATCTAATCCCAGTTTG
# Using the command above the in silico PCR will add the amplicon in the --output just if either, the forward and reverse, primer-binding region is found.
# In my case,  Results | Extracted 22371 amplicons from 883769 sequences (2.53%), just 2.5% of the sequences was added in output 
# Using the parameter --relaxed the amplicon will be added in --output if the forward OR reverse primer-binding region will be finded.
crabs --in-silico-pcr --input crabs_mitofish.txt --output insilico.txt --forward AAACTCGTGCCAGCCACC --reverse GGGTATCTAATCCCAGTTTG --relaxed
# Results | 7582 amplicons were extracted by only the forward or reverse primer (25.31%)

# STEP 5 - RETRIEVE AMPLICONS WITHOUT PRIMER-BINDING REGIONS - this step take some time 
crabs --pairwise-global-alignment --input crabs_mitofish.txt --amplicons insilico.txt --output aligned.txt --forward AAACTCGTGCCAGCCACC --reverse GGGTATCTAATCCCAGTTTG --size-select 10000 --percent-identity 0.95 --coverage 95

# STEP 6 - CURATE AND SUBSET THE LOCAL DATABASE VIA MULTIPLE FILTERING PARAMETERS 
# 6.1 - Dereplicade - Remove identical reference barcodes to speed up taxonomy assignment, as well as improve taxonomy assignment results.
crabs --dereplicate --input aligned.txt --output dereplicated.txt --dereplication-method 'unique_species'

# 6.2 - Filter - six parameters to filter 
# --minimum-length: minimum sequence length for an amplicon to be retained in the database;
# --maximum-length: maximum sequence length for an amplicon to be retained in the database;
# --maximum-n: discard amplicons with N or more ambiguous bases (N);
# --environmental: discard environmental sequences from the database;
# --no-species-id: discard sequences for which no species name is available;
# --rank-na: discard sequences with N or more unspecified taxonomic levels.
crabs --filter --input dereplicated.txt --output filtered.txt --minimum-length 100 --maximum-length 200 --maximum-n 1 --environmental --no-species-id --rank-na 2

# 6.3 - Subset - remove reference barcodes from taxonomic groups not of interest to the research and known erroneous sequences.
crabs --subset --input filtered.txt --output subset.txt --include 'Chordata'

# STEP 7 - EXPORT THE LOCAL DATABASE
# There are 7 formats that you can export the ref database, I'm gonna export in rdp format to use in DADA2
# --export-format 'sintax': The SINTAX classifier is incorporated into USEARCH and VSEARCH;
# --export-format 'rdp': The RDP classifier is a standalone program widely used in microbiome studies;
# --export-format 'qiime-fasta' and --export-format 'qiime-text': Can be used to assign a taxonomic ID in QIIME and QIIME2;
# --export-format 'dada2-species' and --export-format 'dada2-taxonomy': Can be used to assign a taxonomic ID in DADA2;
# --export-format 'idt-fasta' and --export-format 'idt-text': The IDTAXA classifier is a machine learning algorithm incorporated in the DECIPHER R package;
# --export-format 'blast-notax': Creates a local BLAST reference database for blastn and megablast where the output does not provide a taxonomic ID, but lists the accession number;
# --export-format 'blast-tax': Creates a local BLAST reference database for blastn and megablast where the output provides both the taxonomic ID and accession number.
crabs --export --input subset.txt --output 12s_mitofishdb_rdp.fasta --export-format 'rdp'

# STEP 8 - POST-PROCESSING FUNCTIONS TO EXPLORE AND PROVIDE A SUMMARY OVERVIEW OF THE LOCAL REFERENCE DATABASE
# CRABS can run five post-processing functions to explore and provide a summary overview of the local reference database, including 
# (i) --diversity-figure, (ii) --amplicon-length-figure, (iii) --phylogenetic-tree, (iv) --amplification-efficiency-figure, and (v) --completeness-table

crabs --diversity-figure --input subset.txt --output diversity-figure.png --tax-level 4 # (i)
crabs --amplicon-length-figure --input subset.txt --output amplicon-length-figure.png --tax-level 4 # (ii)
crabs --phylogenetic-tree --input subset.txt --output phylo --tax-level 4 --species 'Carcharodon carcharias+Squalus acanthias' # (iii)
crabs --amplification-efficiency-figure --input crabs_mitofish.txt --amplicons subset.txt --forward AAACTCGTGCCAGCCACC --reverse GGGTATCTAATCCCAGTTTG --output amplification-efficiency.png --tax-group Chordata # (iv)
crabs --completeness-table --input subset.txt --output completeness.txt --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp --species 'Carcharodon carcharias+Squalus acanthias' # (v)
