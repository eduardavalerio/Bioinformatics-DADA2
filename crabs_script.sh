##############################################################################
##### Creating Reference databases for Amplicon-Based Sequencing (CRABS) #####
##### Script to create a reference sequences database to Taxonomic ###########
##### Assignment step in DADA2 ###############################################
##### Created by Eduarda V. de Jesus - July 2025 #############################
##############################################################################

# This part is in bash
conda create -n crabs_env python=3.8
conda activate crabs_env
pip install crabs

# Preparing input files 
# Fasta file -> sequences from MitoFish 
# Taxonomy file -> from NCBI or BOLD

# Combine multiple files into a unique one
# Fasta file 
cat *.fa > MitoFish.fasta 
# Taxonomy file 
cat *.gbk > MitoFish.gbk

# From NCBI GenBank Files 
# Extract accessions + taxonomy from .gbk
bioawk -c genbank '{print $name, $qualifiers["organism"]}' MitoFish.gbk | awk '{print $1 "\t" $2}' > taxmap.tsv

# Example Output 
# NC_006919.1    Sundasalanx mesops

# Fetch full lineage (Kingdom→Species)
cut -f2 taxmap.tsv | taxonkit name2taxid | taxonkit lineage | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" > lineages.txt
paste <(cut -f1 taxmap.tsv) lineages.txt > taxmap_final.tsv

# Add full taxonomy ranks using taxonkit
# Fetch full lineage (Kingdom→Species)
cut -f2 taxmap.tsv | taxonkit name2taxid | taxonkit lineage | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" > lineages.txt
paste <(cut -f1 taxmap.tsv) lineages.txt > taxmap_final.tsv

# To ensure the file is tab-separated
head -n1 taxmap.tsv | cat -A  # Should show `^I` for tabs

# Validation
# Count incomplete taxonomies
awk -F'\t' '$2 !~ /;/ {print}' taxmap.tsv | wc -l

# Python
# Run CRABS to generate the .rdp
crabs build_database \
  --input MitoFish.fasta \
  --taxonomy taxmap.tsv \
  --output crabs_db \
  --format rdp \
  --blocks accession,taxonomy

# Filter to Teleo2 primer region 
crabs sequence_extraction \
  --input crabs_db.fasta \
  --primers "ACACCGCCCGTCACTCT" \
  --output teleo2_db.fasta

# Validate output
head -n 1 crabs_db.rdp
