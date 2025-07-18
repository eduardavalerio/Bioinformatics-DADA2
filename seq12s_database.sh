###################################################################################################
####### Script to create a Reference Sequences File to use in Taxonomic Assignment in DADA2 #######
##### Created by Eduarda V. de Jesus - July 2025 ##################################################
################################################################################################### 


# Install references sequencences in MitoFish database 
# We're gonna use the all mitogenomes files 

# Combine multiple fasta files into a unique one
cat *.fa > MitoFish.fasta

# MitoFish sequences often include species names in headers, but you need full taxonomic ranks (Kingdom → Species) in .rdp files 
# Use NCBI Taxonomy to map species name to lineages
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz

# Proceed in python
ipython


from ete3 import NCBITaxa
ncbi = NCBITaxa()

def get_lineage(species_name):
    try:
        taxid = ncbi.get_name_translator([species_name])[species_name][0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        return [names[taxid] for taxid in lineage if ranks[taxid] in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]]
    except:
        return None

# The .rdp format requires:
# >SequenceID;tax=Kingdom;Phylum;Class;Order;Family;Genus;Species
# ACACCGCCCGTC...

import re

# Input: MitoFish FASTA (headers like ">Genus_species")
fasta_file = "MitoFish_12S.fasta"
output_file = "MitoFish_12S.rdp"

with open(fasta_file) as fin, open(output_file, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            # Extract species name (e.g., ">Takifugu_rubripes")
            species = re.sub(r'>| .*', '', line.strip()).replace("_", " ")
            # Get lineage (using ete3 or pre-resolved taxonomy)
            lineage = get_lineage(species)  # e.g., ["Animalia", "Chordata", ...]
            if lineage:
                fout.write(f">{species.replace(' ', '_')};tax={';'.join(lineage)}\n")
            else:
                fout.write(f">{species.replace(' ', '_')};tax=unknown;unknown;unknown;unknown;unknown;unknown;unknown\n")
        else:
            fout.write(line)


# In bash
# Extract accessions (adjust delimiter if needed)
awk -F'|' '/^>/ {print $4}' MitoFish_12S.rdp > accessions.txt

# This step takes a while...
# Install EDirect if needed
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Fetch taxonomy for each accession
while read acc; do
  efetch -db nuccore -id "$acc" -format docsum | \
    xtract -pattern DocumentSummary -element TaxId | \
    xargs -I {} efetch -db taxonomy -id {} -format scientific_name
done < accessions.txt > taxonomy.txt

# Merge accessions with taxonomy
paste -d ";" accessions.txt taxonomy.txt > acc_tax_map.txt

# Replace headers in .rdp file
awk 'BEGIN {FS=";"; OFS=";"} 
     NR==FNR {tax[$1]=$2; next} 
     /^>/ {split($1, a, "|"); acc=a[4]; $0=">" acc ";tax=" tax[acc]} 
     1' acc_tax_map.txt MitoFish_12S.rdp > MitoFish_12S_right.rdp


# Verify Output
head -n 1 MitoFish_12S_right.rdp


# Trim to Teleo2 primer 
cutadapt -g ^ACACCGCCCGTCACTCT... -o MitoFish_Teleo2.rdp MitoFish_12S_right.rdp
# Remove duplicates
cd-hit-est -i MitoFish_Teleo2.rdp -o MitoFish_Teleo2_dedup.rdp -c 1.0

#MitoFish_Teleo2_dedup.rdp is ready to be used in DADA2!
