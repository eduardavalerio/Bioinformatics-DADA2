# Cutadapt pipeline (step before DADA2)
# Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
# https://cutadapt.readthedocs.io/en/stable/index.html

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
#should show the cutadapt version number 
#5.0

cutadapt -a <primer-forward> -o lib1_forward_trimmed.fastq input1.fq.gz
cutadapt -g <primer-reverse> -o lib1_reverse_trimmed.fastq input2.fq.gz


