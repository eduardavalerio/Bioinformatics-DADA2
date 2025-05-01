# Cutadapt pipeline (step before DADA2)
# Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
# https://cutadapt.readthedocs.io/en/stable/index.html

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
# Should show the cutadapt version number 
# 5.0

cutadapt -a AAACTCGTGCCAGCCACC -A GGGTATCTAATCCCAGTTTG -o lib1_1_trimmed.fastq -p lib1_2_trimmed.fastq --minimum-length=130 --quality-cutoff=30 --max-n=0 lib1_1.fq.gz lib1_2.fq.gz

# After this step -> DADA2 R pipeline


