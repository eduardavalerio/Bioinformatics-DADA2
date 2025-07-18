# Install references sequencences in MitoFish database 
# We're gonna use the all mitogenomes files 

# Combine multiple fasta files into a unique one
cat *.fa > MitoFish.fasta

# MitoFish sequences often include species names in headers, but you need full taxonomic ranks (Kingdom → Species) in .rdp files 
# Use NCBI Taxonomy to map species name to lineages
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz

# Proceed in the script seq12s_database.py...
