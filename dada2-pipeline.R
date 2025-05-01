#####################################################
### DADA2 R script - eDNA bioinformatics analysis ###
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

