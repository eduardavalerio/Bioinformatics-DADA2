# Cutadapt pipeline (step before DADA2)
# Demultiplex step 
# https://cutadapt.readthedocs.io/en/stable/index.html
# https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
# cutadapt version 5.0

# To create the barcode.fasta file just follow convert_fasta.py 
cutadapt -e 1 -m 1 -g file:tags_lib1.fasta -o {name}_R1.fastq.gz -p {name}_R2.fastq.gz EV_Lib1_EKDL240035764-1A_22F3HCLT4_L6_1.fq.gz EV_Lib1_EKDL240035764-1A_22F3HCLT4_L6_2.fq.gz

# -e 1 = at least one mismatch in barcode
# -m 1 = remove reads that the lenght is zero. 
# https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
# https://github.com/benjjneb/dada2/issues/159

for file in *_R1.fastq.gz; do echo "$file: $(cat $file | wc -l) reads"; done

# remove fq.gz files with low data
rm -f *.fq.gz


