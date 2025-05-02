# Cutadapt pipeline (step before DADA2)
# Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
# https://cutadapt.readthedocs.io/en/stable/index.html

# Install cutadapt with conda in a new environment 
conda create -n cutadapt cutadapt
conda activate cutadapt 

cutadapt --version
# Should show the cutadapt version number 
# 5.0

# To demultiplex and remove adapters using cutadapt I'm gonna use the NGSfilter, the same as I use in OBITools. 
# Convert NGSfilter for cutadapt accepted file
# Extract SampleID and first barcode (e.g., 'aaaaaaaa' from 'aaaaaaaa:aaaaaaaa')
awk '{print $2 "," $3}' ngsfilter_teleo2_lib1.txt | cut -d ':' -f1 > barcodes_cutadapt.csv
# output barcodes_cutadapt.csv
# 1ALV231.1S,aaaaaaaa
# 2ALV231.2S,aaaaaccc
# 3ALV231.3S,gcaatttt
# 4ALV231.1P,ggacaaac

cutadapt -g ^file:barcodes.csv -G ^GGGTATCTAATCCCAGTTTG --front AAACTCGTGCCAGCCACC --pair-adapters --minimum-length=100 --discard-untrimmed --quality-cutoff=30 --max-n=0 -o lib1_R1.fastq -p lib1_R2.fastq input_R1.fq.gz input_R2.fq.gz


# After this step -> DADA2 R pipeline


