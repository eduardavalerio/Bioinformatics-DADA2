# Cutadapt pipeline (step before DADA2)
# Demultiplex step 
# https://cutadapt.readthedocs.io/en/stable/index.html

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
# Should show the cutadapt version number 
# 5.0

# cutadapt -e 1 -m 1 -g file:barcode.fasta -o output_r1.fastq.gz -p output_r2.fastq.gz input_1.fq.gz input_2.fq.gz
cutadapt -e 1 -m 1 -g file:barcode.fasta -o lib1_R1.fastq.gz -p lib1_R2.fastq.gz EV_Lib2_EKDL240035960-1A_22F7VJLT4_L4_1.fq.gz EV_Lib2_EKDL240035960-1A_22F7VJLT4_L4_2.fq.gz

# -e 1 = at least one mismatch in barcode
# -m 1 = remove reads that the lenght is zero. 
# https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads
# https://github.com/benjjneb/dada2/issues/159

for file in *_R1.fastq.gz; do echo "$file: $(cat $file | wc -l) reads"; done

##remover arquivos com poucos dados fq.gz

rm -f *.fq.gz


