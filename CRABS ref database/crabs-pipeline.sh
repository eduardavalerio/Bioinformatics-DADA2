# Create virtual environment for CRABS and dowload all packages required 
conda create -n crabs-env
conda activate crabs-env
conda install -c conda-forge xopen=1.6.0
conda install -c bioconda cutadapt=4.4
conda install -c bioconda vsearch
pip install crabs 

# Check the installation
crabs -h

# Before start check if you're in the right paste to save all the files 
# I'm gonna use the MitoFish database and NCBI 

# STEP 1 - DOWNLOAD SEQUENCE DATA 
crabs db_download -s mitofish -o mitofish.fasta 

crabs db_download -s ncbi -db nucleotide -q '12S[All Fields] AND ("Chondrichthyes"[Organism] OR "Dipnomorpha"[Organism] 
OR "Actinopterygii"[Organism] OR "Myxini"[Organism] OR "Hyperoartia"[Organism] OR "Coelacanthimorpha"[Organism] OR fish[All Fields])'
-o ncbi12Sfish.fasta -e eduarda.jesus@usp.br

# STEP 2 - MERGE DOWLOADED SEQUENCES
crabs db_merge -o merged.fasta -u yes -i ncbi12Sfish.fasta mitofish.fasta

# STEP 3 - EXTRACT AMPLICONS THROUGH IN SILICO PCR
# In my case, this step I was in troubles bc was returning the error - OverflowError: FASTA/FASTQ record does not fit into buffer
# I filtered sequences longer than 10kb (PCR typically doesn't amplify >10kb)
conda install -c bioconda seqkit
seqkit seq -M 10000 merged.fasta > filtered.fasta 

crabs insilico_pcr -i filtered.fasta -o merged_insilico.fasta -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG #Teleo2 primer fwr and rev

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
