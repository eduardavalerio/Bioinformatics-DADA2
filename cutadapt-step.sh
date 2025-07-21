# Cutadapt pipeline (step before DADA2)
# Demultiplex step 
# https://cutadapt.readthedocs.io/en/stable/index.html

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
# Should show the cutadapt version number 
# 5.0

cutadapt -g CACTCTTTCCCTACACGACGCTCTTCCGATCTTAGARTCTTTGAACGCAAATGGC -G GTGACTGGAGGTTCAGACGTGTGCTCTTCCGATCTGCTTATTAATATGCTTAAATTCA 
-o AI-23-B10_S76_L001_R1_001.trimmed.fastq -p AI-23-B10_S76_L001_R2_001.trimmed.fastq AI-23-B10_S76_L001_R1_001.fastq AI-23-B10_S76_L001_R2_001.fastq --discard-untrimmed  -e 0.12

# -e 1 = at least one mismatch in barcode
# -m 1 = remove reads that the lenght is zero. 
# https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
# https://github.com/benjjneb/dada2/issues/159

for file in *_R1.fastq.gz; do echo "$file: $(cat $file | wc -l) reads"; done

##remover arquivos com poucos dados fq.gz

rm -f *.fq.gz


