#####################################################
### DADA2 R script - eDNA bioinformatics analysis ###
### SCRIPT prepared by Matthieu Leray ###############
### Modified by Eduarda Valerio de Jesus ############
### May 2025 ########################################
#####################################################

# Steps followed

# (1) Assemble metadata #Not done in this manner anymore
# (2) Filter raw reads & Infer Amplicon Sequence Variants (ASVs)
# (3) Assign taxonomy
# (4) Create Phyloseq object
# (5) Remove contaminants, control and outlier samples
#####These alternative extra steps have been removed because they required considerable more time but with a very similar result
# (6) ASV curation ##Only needed if you see diferences in clustering after the vsearch/Lulu step
# (7) Taxonomic assignments with BLAST
# (8) Make a phylogenetic tree
# (9) Create curated Phyloseq object

# Install packages 
install.packages("ggplot2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
if (!require("remotes")) install.packages("remotes")  
remotes::install_github("vmikk/metagMisc")
if (!require("remotes")) install.packages("remotes")
remotes::install_github("tobiasgf/lulu")
install.packages("stringr")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DECIPHER")
install.packages("phangorn")
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mikemc/speedyseq")

# Load packages
library(dada2)
library(ggplot2)
library(phyloseq)
library(decontam)
library(metagMisc)
library(lulu)
library(stringr)
library(DECIPHER)
library(phangorn)
library(speedyseq)

##Creating filepaths to data
#after cutadapt step
path <- "/Volumes/Extreme Pro/uk-data/dada2/lib1/"  #Where the trimmed sequences are
#
#path <- "/Users/quinteroh/Library/CloudStorage/OneDrive-SmithsonianInstitution/eCSI 1st Library/Trimmed_sequences/article_filtered/"
head(list.files(path)) #eventually to check if the path works

##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files(path, pattern = "lib1_1_trimmed.fastq"))
fnRs <- sort(list.files(path, pattern = "lib1_2_trimmed.fastq"))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <-file.path(path, fnFs)
fnRs <-file.path(path, fnRs)

#plotting quality profiles

qprofile_fwd <- plotQualityProfile(fnFs, aggregate = TRUE)
print(qprofile_fwd) + ggtitle("Forward")

qprofile_rev <- plotQualityProfile(fnRs, aggregate = TRUE)
print(qprofile_rev) + ggtitle("Reverse")

#placing filtered files in a new filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Check read length distribution before filtering
library(ShortRead)
hist(width(readFastq(fnFs[1])), breaks = 50, main = "Forward Read Lengths")
hist(width(readFastq(fnRs[1])), breaks = 50, main = "Reverse Read Lengths")

#filtering and trimming, here truncation at 180 (Fwd) and 150 (Rev) bp,  
#2expected errors max (N discarded automatically)  #Truncation length depends on what the qprofile graphs showed
out <- filterAndTrim(
  fnFs, filtFs,          # Input and output files (forward)
  fnRs, filtRs,          # Input and output files (reverse)
  truncLen = c(150, 140), # Forward truncated at 160 bp, reverse at 140 bp -> this changes according to reads length (~150bp)
  maxN = 0,              # Discard reads with ambiguous bases (N)
  maxEE = c(3, 5),       # Allow ≤3 expected errors (forward), ≤5 (reverse)
  truncQ = 10,           # Trim at quality score ≤10 (less aggressive than Q=2)
  minLen = 50,           # Keep reads ≥50 bp after trimming
  rm.phix = TRUE,        # Remove PhiX control reads
  compress = TRUE,       # Output compressed .fastq.gz files
  multithread = TRUE     # Use multiple cores for speed
)

head(out) #eventually to check how filtering and trimming worked

#If not all samples had sequences after filtering, keep only the ones that do
exists0<- file.exists(filtFs) & file.exists(filtRs) #Check which samples were written in the directory
#Keep only the samples with reads after filtration
filtFs <- filtFs[exists0]
filtRs <- filtRs[exists0]

#learning error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
#plotting errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
