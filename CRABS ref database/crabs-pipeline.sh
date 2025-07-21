# Create virtual environment for CRABS and dowload all packages required 
conda create -n crabs-env
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
--email eduarda.jesus@usp.br --database nucleotide
# Taxonomy from NCBI
crabs --download-taxonomy --exclude 'acc2taxid, taxdump'
crabs --download-taxonomy --output ncbi_taxonomy #ncbi_taxonomy is a subdirectory that I crated -> mkdir ncbi_taxonomy 

# STEP 2 - IMPORT DOWNLOADED DATA INTO CRABS FORMAT
# mitofish
crabs --import --import-format mitofish --input mitofish.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp
--acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_mitofish.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'
# ncbi
crabs --import --import-format ncbi --input nci.fasta --names ncbi_taxonomy/names.dmp --nodes ncbi_taxonomy/nodes.dmp
--acc2tax ncbi_taxonomy/nucl_gb.accession2taxid --output crabs_ncbi.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'

# STEP 3 - MERGE DOWLOADED SEQUENCES
# This step is just when sequence data from multiple online repositories are downloaded. After import the data in crabs format step.



















#crabs --download-ncbi --query '("Chondrichthyes"[Organism] OR "Dipnomorpha"[Organism] OR "Actinopterygii"[Organism] 
#OR "Myxini"[Organism] OR "Hyperoartia"[Organism] OR "Coelacanthimorpha"[Organism] OR Fish[All Fields]) AND 12S[All Fields]' 
#--output ncbi_fish12s.fasta --email eduarda.jesus@usp.br --database nucleotide

# STEP 2 - MERGE DOWLOADED SEQUENCES
crabs db_merge -o merged.fasta -u yes -i ncbi12Sfish.fasta mitofish.fasta

# STEP 3 - EXTRACT AMPLICONS THROUGH IN SILICO PCR
# In my case, this step was returning the error - OverflowError: FASTA/FASTQ record does not fit into buffer
# I filtered sequences longer than 10kb (PCR typically doesn't amplify >10kb)
conda install -c bioconda seqkit
seqkit seq -M 10000 merged.fasta > filtered.fasta 

crabs insilico_pcr -i filtered.fasta -o merged_insilico.fasta -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG #Teleo2 primer fwr and rev
# found primers in 7 sequences -> Is this a problem?

# STEP 4 - EXTRACT AMPLICONS THROUGH PAIRWISE GLOBAL ALIGNMENTS
crabs pga -i merged.fasta -o merged_insilicopga.fasta -db merged_insilico.fasta -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG

# STEP 5 - ASSIGNING A TAXONOMIC ID
crabs db_download -s taxonomy
crabs assign_tax -i merged_insilicopga.fasta -o merged_insilicopgatax.tsv -a nucl_gb.accession2taxid -t nodes.dmp -n names.dmp -w yes

head -n 4 merged_insilicopgatax.tsv
# Output
# OR575552	305515	Eukaryota	Chordata	Chondrichthyes	Squaliformes	Dalatiidae	Euprotomicrus	Euprotomicrus_bispinatus	ACTTAAATTAATTATGTAAAATTTTTAACTCTCCGGAGAAAAACCACCCATATAATACCCCTAATTTAACTGTTTTTGGTTGGGGTGACCAGGGGGAAAAAATTATCCCCCCCATCGATTGAGTACTCAGTACTTAAAAATTAGAATGACAACTCTAATTAATAAAACATTTACCGAAACATGACCCAGAATTTATTCTGATCAATGAACCA
# OR582714	862652	Eukaryota	Chordata	Chondrichthyes	Squaliformes	Etmopteridae	Etmopterus	Etmopterus_bigelowi	ACTTAAATTAATCATATAAACTATTAACCCACGGGAATAAATTACATATATACCCCTAATTTAACTGTTTTTGGTTGGGGCGACCAAGGGGGAGAAAAAATCCCCCTCATCGATTGAGTACTTAGTACTTAAAAATTAGAACGACAGTTCTTATTAATGAAATATTTAACGAAAAATGACCCAGTTTTTCTGATCAATGAACCA
# OR582700	263691	Eukaryota	Chordata	Chondrichthyes	Carcharhiniformes	Scyliorhinidae	Apristurus	Apristurus_kampae	ACTTAGACTAATTATGTAATTTTTTTCCGCCTGTGGGTAAAAACAAAAATATAATATTTCTAGTTTAATTGTTTTTGGTTGGGGTGACCAAGGGGAAAAACAAATCCCCCTTATCGACCAAGTACTCAGTACTTAAAAATTAGAGCGACAGCTCTAATCAATAAAACATTTATCGAAAAATGACCCAGGATTTCCTGATCAATGAACCA
# OR582686	671160	Eukaryota	Chordata	Chondrichthyes	Squaliformes	Etmopteridae	Etmopterus	Etmopterus_gracilispinis	ACTTAAATTAATTATGTAAAACTACTAATCCACGGAAATAAACTATTTATATAATATTTCTAATTTAACTGTTTTTGGTTGGGGTGACCGAGGGGAAAAGAAAATCCCCCTCATCGATTGAGTACTTAGTACTTAAAAATTAGAACGACAGTTCTTATTAATAAAATATTTAACGAAAAATGACCCAGTTTTTCTGATCAATGAACCA

# In my case the 'Eukaryota' column is empty -> 'nan'

# STEP 6 - DEREPLICATING THE DB
crabs dereplicate -i merged_insilicopgatax.tsv -o merged_derep.tsv -m uniq_species

# STEP 7 - REF DB CLEANUP 
crabs seq_cleanup -i merged_derep.tsv -o merged_derepcleaned.tsv -e yes -s yes -na 1

# STEP 8 - EXPORT REF DB 
# I'm gonna export ref db in .rdp format
crabs tax_format -i merged_derepcleaned.tsv -o db12S_rdp.fasta -f rdp
head -n 4 db12S_rdp.fasta

# Output
# >OR575552;tax=d:Eukaryota,p:Chordata,c:Chondrichthyes,o:Squaliformes,f:Dalatiidae,g:Euprotomicrus,s:Euprotomicrus_bispinatus
# ACTTAAATTAATTATGTAAAATTTTTAACTCTCCGGAGAAAAACCACCCATATAATACCCCTAATTTAACTGTTTTTGGTTGGGGTGACCAGGGGGAAAAAATTATCCCCCCCATCGATTGAGTA
# >OR582714;tax=d:Eukaryota,p:Chordata,c:Chondrichthyes,o:Squaliformes,f:Etmopteridae,g:Etmopterus,s:Etmopterus_bigelowi
# ACTTAAATTAATCATATAAACTATTAACCCACGGGAATAAATTACATATATACCCCTAATTTAACTGTTTTTGGTTGGGGCGACCAAGGGGGAGAAAAAATCCCCCTCATCGATTGAGTACTTAG

# The output that I had was different than expected 
# >LC878664	root;nan;Chordata;Actinopteri;Cypriniformes;Gobionidae;Hemibarbus;Hemibarbus_labeo
# GCGGTTAAACGAGAGGCCCTAGTTGATATTATTACGGCGTAAAGGGTGGTTAAGGAAAGCACAACAATAAAGCCGAATGGCCCTCTGGCCGTCATACGC

# Is it ok to do that?
awk 'BEGIN{FS="\t|;";OFS=";"} 
     /^>/ {
       print $1 ";tax=d:Eukaryota,p:" $4 ",c:" $5 ",o:" $6 ",f:" $7 ",g:" $8 ",s:" $9
     }
     !/^>/ {
       print
     }' db12S_rdp.fasta > 12S_fish_rdp.fasta

# Output (same as was expected)
# >LC878664;tax=d:Eukaryota,p:Chordata,c:Actinopteri,o:Cypriniformes,f:Gobionidae,g:Hemibarbus,s:Hemibarbus_labeo
# GCGGTTAAACGAGAGGCCCTAGTTGATATTATTACGGCGTAAAGGGTGGTTAAGGAAAGCACAACAATAAAGCCGAATGGCCCTCTGGCCGTCATACGCTTCCAGGTATCCG
