###______________________________________________________________________### 
###                                                                      ###
###                 SAO PAULO STATE UNIVERSITY (UNESP)                   ###
###               MARINE ECOLOGY AND GENOMICS LAB (LAGEM)                ###
###                  COASTAL ISLANDS FROM SP, BRAZIL                     ###
###                       FAPESP (2019/10201-0)                          ###
###         Script optimized by Matthieu Leray and Isabela Buxbaum       ###
###     Minor modifications by Rodrigo R. Domingues and Eduarda Valério  ###
###                                                                      ###
###______________________________________________________________________###



# Sections 
# (1) Set up
# (2) Sequencing Data Processing 
# (3) Inferring Amplification Sequence Variants with DADA 
# (4) Taxonomic Assignment
# (5) Create phyloseq object 
# (6) Remove Contaminants 
# (7) More Filtering 
# (8) Exploratory graphs
# (9) Troubleshooting 
### 1. Minoverlap
### 2. Dada assignment
### 3. Create phyloseq
### 4. Decontam
### 5. OTU clustering 
### 6. OTU curation 
### 7. Insect 
### 8. New phyloseq


####################
#### (1) Set up ####
####################

#install.packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")
BiocManager::install(c("phyloseq", "Biostrings"))
BiocManager::install("dada2")
BiocManager::install("DECIPHER")
install.packages("devtools")
devtools::install_github("tobiasgf/lulu")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
devtools::install_github("mikemc/speedyseq")
devtools::install_github("vmikk/metagMisc")
install.packages("furrr")
install.packages("tidyverse")
install.packages("insect")

#activating and checking packages
library(devtools)
library(dada2)
library(tidyverse)
library(phyloseq)
library(decontam)
library(DECIPHER)
library(speedyseq)
library(metagMisc)
library(lulu)
library(pairwiseAdonis)

# Save your current workspace to a file - just right before close R 
save.image("script_12S.RData")

# Later, load it to restore all variables
load("script_12S.RData")

# set work directory
library(furrr)
plan(multisession, workers = 2)  # two parallel processes, safe for 8GB
options(timeout = 7200)
getwd()
setwd("/Volumes/LAGEM/uk-data/dada2/lib1") # change to your path


########################################
#### (2) Sequencing Data Processing ####
########################################

# creating file path to data
path <- "/Volumes/LAGEM/uk-data/dada2/lib1"
list.files(path)

# check list files, set path, and rename samples
# file preparation

# extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files (path, pattern = "_R1.fastq.gz")) 
fnRs <- sort(list.files (path, pattern = "_R2.fastq.gz"))
sample.names <- sapply(strsplit(fnFs, "_"), '[', 1)
sample.names
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# plotting quality profiles - you can change the parameters to make better the images
tiff("12S_plotquality_forward_sample_aggregate.tiff", width = 6, height = 6, units = "in", res = 300)
qprofile_fwr <- print (plotQualityProfile(fnFs, aggregate = TRUE) 
                       + ggtitle ("Forward"))
dev.off()

tiff("12S_plotquality_reverse_sample_aggregate.tiff", width = 6, height = 6, units = "in", res = 300)
qprofile_rev <- print (plotQualityProfile(fnRs, aggregate = TRUE) 
                       + ggtitle ("Reverse"))
dev.off()

print(qprofile_fwr)
print(qprofile_rev)

# placing filtered files in a new filtered subdirectory
filtFs <- file.path (path, "filtered", paste0 (sample.names, "_F_filt.fastq"))
filtRs <- file.path (path, "filtered", paste0 (sample.names, "_R_filt.fastq"))
names (filtFs) <- sample.names
names (filtRs) <- sample.names

# filtering and trimming, here truncation at 300 (Fwd) and 270 (Rev) bp,
# 2expected errors max (N discarded automatically)
# trimLeft will remove primer + overhang 

# primerForward_Teleo2 <- AAACTCGTGCCAGCCACC - 18 characters
# primerReverse_Teleo2 <- GGGTATCTAATCCCAGTTTG - 20 characters

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(140,150), # inspect your quality profile
                     trimLeft=c(18,20),
                     maxN=0, maxEE=c(2,2), 
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)  # use FALSE in Windows, TRUE in macOS and linux

print(out)

saveRDS(out, "out_12S.rds")
out <- readRDS ("out_12S.rds")

#to just consider files that pass the filtering 
exists <- file.exists (filtFs) & file.exists (filtRs)
filtFs <- filtFs [exists]
filtRs <- filtRs [exists]

# Learning error rates 
errF <- learnErrors (filtFs, multithread = TRUE)
#102067184 total bases in 1116605 reads from 3 samples will be used for learning the error rates.

errR <- learnErrors (filtRs, multithread = TRUE)
#123133443 total bases in 908190 reads from 3 samples will be used for learning the error rates.

# Plotting errors
# https://benjjneb.github.io/dada2/tutorial.html
# The error rates for each possible transition (A→C, A→G, …) are shown. Points 
# are the observed error rates for each consensus quality score. The black line 
# shows the estimated error rates after convergence of the machine-learning 
# algorithm. The red line shows the error rates expected under the nominal 
# definition of the Q-score. Here the estimated error rates (black line) are a 
# good fit to the observed rates (points), and the error rates drop with increased
# quality as expected. Everything looks reasonable and we proceed with confidence.

# Some error rates are 0 (especially at high quality scores), and log10(0) = -∞ (infinite)
# These points get automatically filtered out, leading to the warning.

tiff("12S_errF.tif", width = 10, height = 10, units = "in", res = 300)
plotErrors(errF, nominalQ=TRUE) 
dev.off() # Removed 82 rows containing missing values or values outside the scale range

tiff("12S_errR.tif", width = 10, height = 10, units = "in", res = 300)
plotErrors (errR, nominal = TRUE)
dev.off() # Removed 82 rows containing missing values or values outside the scale range

# Dereplicating reads
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
derepFs <- derepFastq (filtFs)
names (derepFs) <- sam.names
derepRs <- derepFastq (filtRs)
names (derepRs) <- sam.names

# Save copies of dereplicated objects 
saveRDS(derepFs, file ="/Volumes/LAGEM/uk-data/dada2/lib1/derepFs.rds")
saveRDS(derepRs, file ="/Volumes/LAGEM/uk-data/dada2/lib1/derepRs.rds")

derepFs <- readRDS ("derepFs.rds")
derepRs <- readRDS ("derepRs.rds")


################################################################
#### (3) Infering Amplification Sequence Variants with dada ####
################################################################

dadaFs <- dada(derepFs, 
               err = errF, 
               multithread=TRUE,
               pool=F) 

dadaRs <- dada(derepRs, 
               err=errR,
               multithread=TRUE,
               pool=F)

dadaFs[[1]]
dadaRs[[1]]

# Save copies of dada objects
saveRDS(dadaFs, file ="/Volumes/LAGEM/uk-data/dada2/lib1/dadaFs_12s.rds")
saveRDS(dadaRs, file ="/Volumes/LAGEM/uk-data/dada2/lib1/dadaRs_12s.rds")

dadaFs <- readRDS ("dadaFs_12s.rds")
dadaRs <- readRDS ("dadaRs_12s.rds") 

# Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
head(mergers[[1]])

# Construct sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#My result:
# [1]    77 3871
# This means in 77 samples we have 3871 different ASVs - is that low?

# Inspect distribution of sequence lenghts 
table (nchar(getSequences(seqtab)))

# Saving files to use in the next part of the workflow
saveRDS(seqtab, file ="/Volumes/LAGEM/uk-data/dada2/lib1/seqtab_12S.rds")

# Save copies of seqtab (sequence table)
write.csv(seqtab, file ="/Volumes/LAGEM/uk-data/dada2/lib1/seqtab+Chim.csv")

# Identifying and removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                     method="consensus", 
                                     multithread=TRUE)

dim (seqtab.nochim) # 77 2799
sum(seqtab.nochim)/sum(seqtab) # 0.9940098

# Saving table without chimeras and downloading as .csv file
saveRDS(seqtab, file = "/Volumes/LAGEM/uk-data/dada2/lib1/seqtab_nochim_12S.rds")
write.csv (seqtab.nochim, file =  "/Volumes/LAGEM/uk-data/dada2/lib1/nochim_seq_12S.csv")
write.table(seqtab.nochim, file = "/Volumes/LAGEM/uk-data/dada2/lib1/nochim_seq_12S.txt", sep = "\t", quote = FALSE, col.names = NA)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                  sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sam.names
head(track)

# This is a great place to do a last sanity check. Outside of filtering, 
# there should no step in which a majority of reads are lost. If a majority 
# of reads failed to merge, you may need to revisit the truncLen parameter 
# used in the filtering step and make sure that the truncated reads span your 
# amplicon. If a majority of reads were removed as chimeric, you may need to 
# revisit the removal of primers, as the ambiguous nucleotides in unremoved 
# primers interfere with chimera identification.

# Saving table
write.csv(track, file = "/Volumes/LAGEM/uk-data/dada2/lib1/tracked_reads_12S.csv")
saveRDS(seqtab.nochim, file = "/Volumes/LAGEM/uk-data/dada2/lib1/seqtab.nochim.rds")


##################################
#### (4) Taxonomic Assignment ####
##################################

# This step I ran using a custom database created using CRABS
# follow: https://github.com/eduardavalerio/Bioinformatics-DADA2/tree/main/CRABS%20ref%20database

nochim.tax.merged <- assignTaxonomy(seqtab.nochim, "/Volumes/LAGEM/uk-data/dada2/12s_merged_rdp.fasta", multithread = TRUE, tryRC = TRUE) #12s_mitofish_rdp.fasta - Crabs database
saveRDS(nochim.tax.merged, file = "/Volumes/LAGEM/uk-data/dada2/lib1/nochim.tax.rds")

taxa.print <- nochim.tax # removing sequences rownames for display only 
rownames(taxa.print) <- NULL
head(taxa.print)

tax_dada = tax_table(nochim.tax)
write.csv(tax_dada, file = "/Volumes/LAGEM/uk-data/dada2/lib1/nochim.tax.csv")

otu_dada = otu_table(seqtab.nochim, taxa_are_rows = FALSE)
write.csv(otu_dada, file = "/Volumes/LAGEM/uk-data/dada2/lib1/seqtab.nochim.csv")

dim(otu_dada)
dim(tax_dada)

# Proximos passos: criar banco de dados mitofish + ncbi e rodar assignTaxonomy dnv usando multithread = TRUE -> antes tava = 8


#________________________________________________________________________PAREI AQUI
####################################
#### (5) Create phyloseq Object ####
####################################

#Creating a phyloseq object - charge necessary libraries
library(phyloseq)
library(tibble)
library(dplyr)
library(phyloseq)
library(Biostrings)

#At the end of this script, I found trouble concluding the taxonomic assignment because the original script did not include the csv with the reference sequences in the phyloseq obj.
#So I added a few more steps in order to create another phyloseq containing the ref seqs, which is the one we will be using afterwards to continue with the analysis

#Additionally, the csv files automatically generated previously by dada2 (ASVtable - nochim_seq_new_names.csv) and insect (TAXtable - taxonomic_assignment) are generated with the wrong format.
#So make sure to create new csv files, containing the same information, naming the necessary columns with 'SampleID', 'ASVS' and excluding any extra (and unnecessary) column/line that you may come across.

#Also, create a new csv file containing all the metadata from your data (it's the SAMtable)

#Creating phyloseq object - reading the 3 csv tables: 
#ASVtable from dada2
#TAXtable from insect
#SAMtable - a new table with the variable information of the samples
asvtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\nochim_seq_new_names_useful.csv")
taxtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\taxonomic_assignment_useful.csv")
samtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\SAMtable_plato2.csv")

#Define rownames
asvtable <- asvtable %>% tibble::column_to_rownames("SampleID")
taxtable <- taxtable %>% tibble::column_to_rownames("ASVS")
samtable <- samtable %>% tibble::column_to_rownames("SampleID")

#Converte all columns to numeric
asvtable[] <- lapply(asvtable, function(x) as.numeric(as.character(x)))

#transform in matrices
asvtable <- as.matrix(asvtable)
taxtable <- as.matrix(taxtable)

#CREATE the phyloseq object
ROHR42.1_obj <- phyloseq(otu_table(asvtable, taxa_are_rows = FALSE), 
                         sample_data(samtable), tax_table(taxtable))

#removings ASVs that are not present in any sample (if any)
ROHR42.1_obj <- prune_taxa(taxa_sums(ROHR42.1_obj) > 0, ROHR42.1_obj)

#SAVE phyloseq object
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\ROHR42.1_obj.rds")

#Verifying what's inside this phyloseq obj that was just created
ROHR42.1_obj
otu_table(ROHR42.1_obj)                  #returns the abundance matrix
otu_table(ROHR42.1_obj)[1:5, 1:5]        #returns the first 5 ASVs x 5 samples

tax_table(ROHR42.1_obj)                  #returns the taxonomic matrix
tax_table(ROHR42.1_obj)[1:5, ]           #returns the first 5 ASVs with taxonomic information

sample_data(ROHR42.1_obj)                #returns the metadata
as.data.frame(sample_data(ROHR42.1_obj)) #returns the metadataas a dataframe

sample_names(ROHR42.1_obj)               #returns the names of the samples
taxa_names(ROHR42.1_obj)                 #returns the names of the ASVs
rank_names(ROHR42.1_obj)                 #returns the taxonomic levels (e.g. Kingdom, Genus)

str(ROHR42.1_obj)
slotNames(ROHR42.1_obj)                  #returns what are the slots present in the phyloseq obj

ROHR42.1_obj@otu_table


#Charge CSV with the reference sequences
seqfile <- "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\reference_plato2.csv"
seqdata <- read.csv(seqfile)

#Verify the header to make sure all columns are compatible
head(seqdata)

#Create a DNAStringSet object with the reference sequences from the csv
asv_seqs <- DNAStringSet(seqdata$Sequence)

#Name sequences according to the ASVs
names(asv_seqs) <- seqdata$ASVS

#Verify if the names of the ASVs match with the taxa_names of the phyloseq obj
all(names(asv_seqs) %in% taxa_names(ROHR42.1_obj)) #Result should be TRUE

#Add sequences to the refseq slot of the phyloseq obj
ROHR42.1_obj@refseq <- asv_seqs

#Verify if the atribution worked
ROHR42.1_obj@refseq[1:5]  # Verifica as primeiras 5 sequências

#Save the new phyloseq obj with the new refseq slot
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\full_ROHR42.1_obj.rds")



#________________________________________________________________________
#
#### Remove contaminants ####
#________________________________________________________________________

#Plot Controls vs true samples

library(decontam)
packageVersion("decontam")

df <- as.data.frame(sample_data(ROHR42.1_obj))
df$LibrarySize <- sample_sums(ROHR42.1_obj)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#Plot: controls (is.neg + TRUE) vs true samples
ggplot(data = df, aes(x=Index, y=LibrarySize, color=is.neg)) + geom_point()  #antes estava color=Sample, mudei para color=is.neg

#Identify contaminants - Prevalence
contamdf.prev <- isContaminant(ROHR42.1_obj, method="prevalence", neg="is.neg")
contaminants <- table(contamdf.prev$contaminant) #identify how many contaminants were detected
write.csv(contamdf.prev, 
          file = "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_contamdf.prev.csv")

#Make phyloseq object of presence-absence in negative controls and true samples
ROHR42.1_obj.pa <- transform_sample_counts(ROHR42.1_obj, function(abund) 1*(abund>0))
ROHR42.1_obj.pa.neg <- prune_samples(sample_data(ROHR42.1_obj.pa)$is.neg == TRUE, ROHR42.1_obj.pa)
ROHR42.1_obj.pa.pos <- prune_samples(sample_data(ROHR42.1_obj.pa)$is.neg == FALSE, ROHR42.1_obj.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ROHR42.1_obj.pa.pos), pa.neg=taxa_sums(ROHR42.1_obj.pa.neg), 
                    contaminants=contamdf.prev$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminants)) + geom_point() + 
  xlab("Prevalence(Negative Controls)") + ylab ("Prevalence(True Samples)")

# The original script attempted to create a new csv file using phyloseq obj, but this is not possible so I changed the code a little bit.
# Now, I'm creating a new csv using a dataframe

df.pa$ASV <- rownames(df.pa)
write.csv(df.pa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\df.pa.csv", row.names = FALSE)

#Save work to this point
save(ROHR42.1_obj, contamdf.prev, file = "my_important_objects.RData")
save.image("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\phyloseq-decontam.RData")


# Create a new phyloseq object without contaminants identified in contamdf.prev
# The original script removed SEQUENCES that were KNOWN as contaminants, since we don't know that explicitly up to this point, I did it using the information we've got so far


contamdf.prev <- isContaminant(ROHR42.1_obj, method="prevalence", neg="is.neg")

#Extract the names od the ASVs determined as contaminants

badTaxa <- rownames(contamdf.prev)[which(contamdf.prev$contaminant)]
length(badTaxa)         #shows how many ASVs are considered contaminantes
head(badTaxa, n = 5)    #shows contaminants

#Create new phyloseq without the contaminants
goodTaxa <- setdiff(taxa_names(ROHR42.1_obj), badTaxa)
ROHR42.1_obj <- prune_taxa(goodTaxa, ROHR42.1_obj)

#Now remove all control samples 
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR4")
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR5")
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR6")

#Look at distribution of library size again
df <- as.data.frame(sample_data(ROHR42.1_obj))
df$LibrarySize <- sample_sums(ROHR42.1_obj)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

##Distribution of library size
ggplot(data = df, aes(x = Index, y = LibrarySize)) + 
  geom_point()

##Distribution of library size per site
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Site)) + 
  geom_point()

##Distribution of library sizer per depth
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Depth)) + 
  geom_point()

##Distribution of library size per season
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Season)) + 
  geom_point()

#### More Filtering ####

#removings ASVs that are not present in any sample (if any)
ROHR42.1_obj <- prune_taxa(taxa_sums(ROHR42.1_obj) > 0, ROHR42.1_obj)

#removing ASVs assigned to Eukaryotes and unassigned (NA) at the kingdom level
ROHR42.1_obj <- subset_taxa(ROHR42.1_obj, kingdom != "Eukaryota")

#removing ASVs of Chlorophyta as it can interfere with the coral analysis, 
#since we expect a high level of upwelling microorganisms on the water
ROHR42.1_obj <- subset_taxa(ROHR42.1_obj, phylum != "Chlorophyta")

#Last clean 
ROHR42.1_obj <- prune_taxa(taxa_sums(ROHR42.1_obj) > 0, ROHR42.1_obj)
ROHR42.1_obj 


#My results:
##phyloseq-class experiment-level object
##otu_table()   OTU Table:          [ 1866 taxa and 93 samples ]:
##sample_data() Sample Data:        [ 93 samples by 22 sample variables ]:
##tax_table()   Taxonomy Table:     [ 1866 taxa by 11 taxonomic ranks ]:
##refseq()      DNAStringSet:       [ 1866 reference sequences ]
##taxa are columns


#Save phyloseq object
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

#Reading the phyloseq object
ROHR42.1 <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

#Checking what's inside the phyloseq obj
ROHR42.1_obj
otu_table(ROHR42.1_obj)                   #shows abundance matrix
otu_table(ROHR42.1_obj)[1:5, 1:5]         #first 5 ASVs x 5 samples

tax_table(ROHR42.1_obj)                   #shows taxonomy matrix
tax_table(ROHR42.1_obj)[1:5, ]            #first 5 ASVs with taxonomic information

sample_data(ROHR42.1_obj)                 #shows metadata
as.data.frame(sample_data(ROHR42.1_obj))  #metadata in dataframe format

sample_names(ROHR42.1_obj)                #sample names
taxa_names(ROHR42.1_obj)                  #ASVs names
rank_names(ROHR42.1_obj)                  #taxonomy levels

str(ROHR42.1_obj)
slotNames(ROHR42.1_obj)

ROHR42.1_obj@refseq

#Check first sequences of the refseq
head(refseq(ROHR42.1_obj))

#Check if the refseq slot contains the sequences before saving
head(refseq(ROHR42.1_obj))

if (is.null(refseq(ROHR42.1_obj))) {
  print("Refseq is NULL. Please add sequences before saving.")
} else {
  print("Refseq contains sequences.")
}

saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

loaded_obj <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")
head(refseq(loaded_obj))

loaded_obj <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_2\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")
loaded_obj

#________________________________________________________________________
#
#### Exploratory Graphs ####
#________________________________________________________________________


#Merge of species that have the same taxonomy at a certain taxonomic rank. 
rh1 <- tax_glom(ROHR42.1_obj, taxrank = "order")
head(tax_table(rh1))

#Removing NA
#Install metagMisc (you will need to install vegan package first as a dependencie of metagMisc)
devtools::install_github("vmikk/metagMisc")
library(metagMisc)
na.rh1 <- phyloseq_rm_na_tax(rh1)

#Draw bar plots
plot_bar(na.rh1, fill = "order")
plot_bar(na.rh1, fill = "order") + geom_bar(aes(color = order, fill = order), 
                                            stat = "identity", position = "stack")
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Location)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Depth)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = Site~Depth)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Season)

#Prepare for ordination plots and do plots
na.rh1 <- prune_samples(sample_sums(na.rh1) > 0, na.rh1)
rh1.ord <- ordinate(na.rh1, "NMDS", "bray")
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Site", 
                title = "Samples by Locality") + geom_point(size = 2)
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Season", 
                title = "Samples by Season") + geom_point(size = 2)
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Depth", 
                title = "Samples by profundidade") + geom_point(size = 2)
rh1_ord2 <- ordinate(na.rh1, "PCoA", "bray")
plot_ordination(na.rh1, rh1_ord2, type = "samples", color = "Site", 
                title = "PCoA for Samples by Locality") + geom_point(size = 2)

#Taxa after filtering of NA results
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 17 taxa and 99 samples ]
#sample_data() Sample Data:       [ 99 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 17 taxa by 5 taxonomic ranks ]

#normalizing number of reads for relative abundance
total = median(sample_sums(na.rh1))
standf = function(x, t=total) round(t * (x/sum(x)))
na.rh1 = transform_sample_counts(na.rh1, standf)
na.rh1

#drawing plots with normalized data
plot_bar(na.rh1, fill = "order")

#merging samples by locality for ploting
sample_data(na.rh1)$Site <- factor(sample_data(na.rh1)$Site)
sum(is.na(sample_data(na.rh1)$Site))
na.rh1 <- subset_samples(na.rh1, !is.na(Site))
na.rh1_locality <- merge_samples(na.rh1, "Site")
plot_bar(na.rh1_locality, fill = "order") 

#plotting heatmap for abundance by locality
otu_table(na.rh1) <- otu_table(otu_table(na.rh1) + 1, taxa_are_rows = TRUE)
(p <- plot_heatmap(na.rh1, "NMDS", "bray", "Site", "order"))

p <- plot_heatmap(na.rh1, method = "NMDS", distance = "bray", sample.label = "Site", taxa.label = "order")
print(p)  # Força a exibição


#Filtering data to consider taxa present at least at 20%
rh1_abund <- filter_taxa(na.rh1, function(x) sum(x > total*0.20) > 0, TRUE)
rh1_abund

#run heatmap for filtered data
(q <- plot_heatmap(rh1_abund, "NMDS", "bray", "Site", "order"))
q <- plot_heatmap(na.rh1, method = "NMDS", distance = "bray", sample.label = "Site", taxa.label = "order")
print(q)

#ordination plots with normalized data
#generating object for nmds with bray curtis and ploting
normalize.ord <- ordinate(na.rh1, "NMDS", "bray")
plot_ordination(na.rh1, normalize.ord, type = "samples", color = "Site", 
                title = "Samples by Locality") + geom_point(size = 3)

#generating object for PCoA with bray curtis and ploting
normal.pcoa <- ordinate(na.rh1, method = "PCoA", distance = "bray")
plot_ordination(na.rh1, normal.pcoa, type = "samples", color = "Site", 
                title = "PCoA for Samples by Locality") +
  geom_point(size = 4)


#phyloseq graphs...



####__________________________________________________________________________
####
#### TROUBLESHOOTING - Repeating the taxonomic assignment ####
####__________________________________________________________________________

## To answer the question of why we don´t have Pocillopora sequences,
## let's extract the sequences obtained by each step.  
## Save derep and dada objects: lines 84, 85, 105 and 106 in this script. 

# Set working directory and read objects
setwd("C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela/Fastq_ROHR42/Alcatrazes_2023/Alcatrazes_2023_raw/filtered/Dereplicated")

derepFs <- readRDS("derepFs.rds")
derepRs <- readRDS("derepRs.rds")
dadaFs <- readRDS("dadaFs_ITS2.rds")
dadaRs <- readRDS("dadaRs_ITS2.rds")

# --> Take a look before merging 
derep_tab <- makeSequenceTable(derepFs, derepRs)

# --> Then, take a look before removing chimeras 
# by recreating or checking csv file saved in line 136 of this script. 


### Trouble shooting ###

## After the first analysis of this data set, none Pocillopora sequences were obtained. 
## Pocillopora genus is the most common taxa composing the Pacific Panamanian reefs assemblage. 
## To add this key taxa to the data set, editions are being done to the script: 

### 1. During the mergers step, the function "minOverlap" is changed from 20 (default) to 1 ####
#  --> This to consider long sequences as the its2 sequence from Pocillopora 
#  --> which can be found in GenBank in a range from 414 - 706. 

merge1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 1, verbose = T)
seqtab.mer1 <- makeSequenceTable (merge1)
dim(seqtab.mer1)

# --> 117 samples --> 45628 asv

table (nchar(getSequences(seqtab.mer1)))

# 220   221   222   224   225   226   227   228   229   230   231   232   233   234   235   236   237   238   239   240   241 
# 49     1     1     6     9    48    12    10     8    10     5    20     8     8     3     9     8    10    18     6    17 
# 242   243   244   245   246   247   248   249   250   251   252   253   254   255   256   257   258   259   260   261   262 
# 21     7    17     4    23    53    15   106    20    29    30    14    17    22    22    28    25    43    32    67   363 
# 263   264   265   266   267   268   269   270   271   272   273   274   275   276   277   278   279   280   281   282   283 
# 176   106   266    54    20    51    44    57   398   200    32    18    13    30    21    33    24    16    43   210    26 
# 284   285   286   287   288   289   290   291   292   293   294   295   296   297   298   299   300   301   302   303   304 
# 32    36    81    46   619     5    37    32    19    17     8    15    27   166    25    31   427     9     3     2     5 
# 306   307   308   309   310   311   312   313   314   315   316   317   318   319   320   321   322   323   324   325   326 
# 182    22    24   344    26    11    23    46     8    13     6    59     1    10    16     5    66    13     2    10    13 
# 327   328   329   330   331   332   333   334   335   336   337   338   339   340   341   342   343   344   345   346   347 
# 40    13     4    20     7     5     9     9     6    28    24    86     7    47    18    18    11    33   116    52    57 
# 348   349   350   351   352   353   354   355   356   357   358   359   360   361   362   363   364   365   366   367   368 
# 24    46   141    93    44    16    38   129    23    30    64     4    12    19    16     4   194    13     4    19    90 
# 369   370   371   372   373   374   375   376   377   378   379   380   381   382   383   384   385   386   387   388   389 
# 15    23     5    21    29    64    43    48   132    52    46     8    37    15     3    18    16    38    24    32     4 
# 390   391   392   393   394   395   396   397   398   399 
#   9    12    18     4    26    52   176  2622  5018 29781


### --> NOTA: para que excel no haga ¡CAPUT! transponemos la tabla

tseqtab.mer1 <- t(seqtab.mer1)

# Save obtained sequences
write.csv (tseqtab.mer1, file = "C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela/Fastq_ROHR42/Alcatrazes_2023/Alcatrazes_2023_raw/filtered/Dereplicated/corrections/seqtab.mer1.csv")
saveRDS(seqtab.mer1, file ="C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela/Fastq_ROHR42/Alcatrazes_2023/Alcatrazes_2023_raw/filtered/Dereplicated/corrections/seqtab.mer1.rds")

# Let's remove the chimeras in it: 

seqtab.mer1.nochim <- removeBimeraDenovo(seqtab.mer1, 
                                         method="consensus", 
                                         multithread=8, 
                                         verbose = T)
dim (seqtab.mer1.nochim)

# Resume:
# --> 117 samples --> 18,124 ASV
# --> 27,504 Chimeric ASV's removed

tseqtab.mer1.nochim <- t(seqtab.mer1.nochim)

write.csv (tseqtab.mer1.nochim, file = 
             "C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela/Fastq_ROHR42/Alcatrazes_2023/Alcatrazes_2023_raw/filtered/Dereplicated/corrections/seqtab.mer1.nochim.csv")
saveRDS(seqtab.mer1.nochim, file =
          "C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela/Fastq_ROHR42/Alcatrazes_2023/Alcatrazes_2023_raw/filtered/Dereplicated/corrections/seqtab.mer1.nochim.rds")


### 2. A different pipeline is used to generate the taxonomy table ####
### --> for creating the phyloseq object

### In this case, run a first assignment with dada2 just to create
### an empty tax table that will be used to run phyloseq and, then, decontam. 
### For this end, I decided to use a "freshly created" TEP Pocillopora database obtained from GenBank with magicBlast

# assign taxonomy to the new no-chimera seq tab 

nochim.tax <- assignTaxonomy(seqtab.mer1.nochim,
                             "C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/TEP_poci.fasta",
                             multithread = 8)
saveRDS(nochim.tax, "C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/merge1.nochim.tax.rds")


taxa.print <- nochim.tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_dada  = tax_table(nochim.tax)
write.csv(tax_dada, file="C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/nochim.tax.csv")

otu_dada  = otu_table(seqtab.mer1.nochim, taxa_are_rows = FALSE)
write.csv(otu_dada, file="C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/seqtab.mer1.nochim.csv")

dim(otu_dada)
dim(tax_dada)

### 3. Create phyloseq ####

meta <- read.table("metadata.csv", header = T, sep = ",", row.names = 1)

ITS2 <- phyloseq(otu_dada, sample_data(meta), tax_dada)
ITS2

saveRDS(ITS2, "C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/ITS2.rds")


### 4. Remove contaminants, controls and outliers with decontam ####

ITS2 <- readRDS("ITS2.rds")
variable.names(ITS2@sam_data)

sample_data(ITS2)$is.neg <- sample_data(ITS2)$sam_type == "Mock"
contamdf.prev2 <- isContaminant(ITS2, method = "prevalence", neg = "is.neg", threshold = 0.5)
contaminants <- table(contamdf.prev2$contaminant)
write.csv(contamdf.prev2, file = "~/2023/Uach/TopicosAvanzados/pipeline/contaminantes.csv")

# Checking contaminants distribution

its2.pa <- transform_sample_counts(ITS2, function(abund) 1*(abund>0))           # presence - abscence 
its2.pa.neg <- prune_samples(sample_data(its2.pa)$sam_type == "Mock", its2.pa)  # extract controls
its2.pa.pos <- prune_samples(sample_data(its2.pa)$sam_type == "eDNA", its2.pa)  # extract samples 
df.pa <- data.frame(pa.pos = taxa_sums(its2.pa.pos), pa.neg = taxa_sums(its2.pa.neg), 
                    contaminant = contamdf.prev2$contaminant)                   # create data frame

# Plotting

ggplot(data = df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Establish contaminant sequences

badTaxa2 = c("GCCCCCTGGTTCTCCAGGGAGCATGTCTGTCTGAGCGTCATGTAACACTAGTACACGCACACTCTTCATGTGTGTGCGGCGTTGAGGCTTCATGGCTTGCTCTTAGTGAGCGAGCTGCGTGCCTTGAAGTGCAGAATATGACTGGGACCCTTGCACTTGCACACGAGTGTGATAGGGGGCAGCAAAAAGAGTAATCTTGTCTGCTTGCCTGTCTCCACGTGCGAGTGCTCGTCGTCTCTTTGGCTTATGCCAAACACTGACTTTGACCTCAGATCAGGCAAGACTACC", 
             "GCTTCTGGTTCCGCCAGGAGCACGTCTGTCTGAGAGTTTACTTTACTGTGCACTGCCCACCCGATGGGTGAGTGTGTTGTGGGGTGTTGTCTGTAGGCCATCGTGTCACGGGCATCCCTCGAAAGGATTGCATCCCTCCCGATTGTGGGAGCCAGCACACAGATGCTGTCGCTGGCTGTCCAGTCGATGGCTGGGGTACTCTGACCTTAACTATCCCTTTAACTGAACCTCAGATCAGGCGAGACTACC", 
             "GCTCTTGGGTTCTCCCAGGAGCATGTCTGTCTGAGTGTCGGATTTTATCGAACGCACTCGTGAGTGCGGTATTGAGGTGTCACGGCGTGTTAGCATTAAATCGCCGTGTCCCTTGAAGGACAGCAGCATTGGACTCGCATTCTCTATTGAAACTTGAGTGCGAAGGCTAAAGATAATATTTTCCTTCGCCCAAACATGCAATAGCAGAAACAAAATCCAATCAAAGCGGTGACAGACAATGACAACATCAAATCTTGACCTCAGATCAGGCGAGATCACC", 
             "GCTCTCGGTAAGGCCGGGAGCACGTCTGCCTGAGTGTCTGATTAAGTTGGAATGTATATTTGGGAGTCGGTTTGGCTTATTCTGCCTTTAGCGTCTCCTGCAATGTGAATAACATTGGCTGGGTTTTTTGCCTTGGTGGCACATTTAGAAGTGAAATCAAGCAAAATGCCCAGGTAAATAATATGGGCCAAGGTGTGCAAGATTTTGATGGACCTCAGATCAGGCGAGAACACC")

goodTaxa2 <- setdiff(taxa_names(ITS2), badTaxa2)
ITS2 <- prune_taxa(goodTaxa2, ITS2)

# Remove control samples 

ITS2 = subset_samples(ITS2, sample_names(ITS2) !="cntrl1-ITS2-PCRI")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="cntrl2-ITS2-PCR2")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="cntrl2-ITS2-PCRI")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="cntrl3-ITS2-PCR2")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="cntrl3-ITS2-PCRI")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3621-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3634-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3639-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3652-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3673-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3882-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3895-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3896-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3913-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML3992-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML4001-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML4167-ITS2-EXT")
ITS2 = subset_samples(ITS2, sample_names(ITS2) !="ML4218-ITS2-EXT")

# Remove any sample that has 0 reads.
ITS2 = subset_samples(ITS2, sample_names(ITS2) != "...")

# removings ASVs that are not present in any sample (if any)
ITS2 <- prune_taxa(taxa_sums(ITS2) > 0, ITS2)
ITS2

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18120 taxa and 99 samples ]
# sample_data() Sample Data:       [ 99 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 18120 taxa by 7 taxonomic ranks ]

saveRDS(ITS2, file = "C:/Users/Viviane/Documents/2023/Uach/TopicosAvanzados/pipeline/its2-free.rds")


# Let's check the distribution of the library size 

df <- as.data.frame(sample_data(ITS2)) 
df$LibrarySize <- sample_sums(ITS2)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(df, aes(x = Index, y = LibrarySize, color= sam_type)) + geom_point() + theme_bw()


### 5. Clusterinng ASVs into 97% OTUs ####

setwd("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU")

nproc <- 8
its2 <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\ROHR42.1_final_obj.rds")
sequences <- taxa_names(its2)
sample_names <- sample_names(its2)

# Verifica se o refseq está presente
if (is.null(refseq(its2))) {
  stop("O slot 'refseq' está vazio. Verifique se o objeto salvo contém as sequências.")
}

# Install decipher 
# nao fiz: dna <- Biostrings::DNAStringSet(sequences)

# fiz: Pegar as sequências diretamente do slot refseq
dna <- refseq(its2)  # <- ESTA É A CORREÇÃO
set.seed(123) # initialize the random number generator
clusters <- DECIPHER::Clusterize(dna, method = "overlap", 
                                 cutoff = 0.03, #97% 
                                 penalizeGapLetterMatches = NA, 
                                 includeTerminalGaps = TRUE, 
                                 processors = nproc)
set.seed(NULL)


# penalize gap-to-letter mismatches once per insertion or deletion, 
# which treats runs of gaps (i.e., indels) as equivalent to a single mismatch
# the calculation of distance will use the entire (global) alignment

# Note: function "merge_taxa_vec" was in the package mikemc/speedyseq 
# NOTE: 99% OTUS gave 6004 ASV's | 97% OTUS gave 2408 ASV's --> I'm choosing 97% otus for the analysis. 
# Given that diversity obtained through both approaches is basically the same. 

# Merge clusters with taxonomy
its2.otu <- merge_taxa_vec(
  its2,
  group = clusters$cluster,
  tax_adjust = 0)

# Checking objects 

its2
#Viviane
# otu_table()   OTU Table:          [ 18120 taxa and 99 samples ]:
# sample_data() Sample Data:        [ 99 samples by 16 sample variables ]:
# tax_table()   Taxonomy Table:     [ 18120 taxa by 7 taxonomic ranks ]:

#Isabela
#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 1430 taxa and 72 samples ]:
#  sample_data() Sample Data:        [ 72 samples by 4 sample variables ]:
#  tax_table()   Taxonomy Table:     [ 1430 taxa by 11 taxonomic ranks ]:
#  refseq()      DNAStringSet:       [ 1430 reference sequences ]
#  taxa are columns


its2.otu #99%

#Viviane
# otu_table()   OTU Table:          [ 6004 taxa and 99 samples ]:
# sample_data() Sample Data:        [ 99 samples by 16 sample variables ]:
# tax_table()   Taxonomy Table:     [ 6004 taxa by 7 taxonomic ranks ]:

#Isabela
#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 608 taxa and 72 samples ]:
#sample_data() Sample Data:        [ 72 samples by 4 sample variables ]:
#tax_table()   Taxonomy Table:     [ 608 taxa by 11 taxonomic ranks ]:
#refseq()      DNAStringSet:       [ 608 reference sequences ]
#taxa are columns

# Take a look at files
tax_dada_tax  = tax_table(its2.otu)
write.csv(tax_dada_tax, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\otu_tax.csv")

tax_dada_otu  = otu_table(its2.otu, taxa_are_rows = FALSE)
write.csv(tax_dada_otu, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\otu_matrix.csv")

##Passo extra para gerar um csv contendo as sequencias de referencia do phyloseq
# 1. Verifique se há uma referência de sequência (refseq)
refseqs_clustered <- refseq(its2.otu)

# 2. Crie um data frame com os nomes dos ASVs (no caso, OTUs) e as sequências
asv_seq_df <- data.frame(
  ASV = names(refseqs_clustered),
  Sequence = as.character(refseqs_clustered)
)

# 3. Escreva para um arquivo CSV
write.csv(asv_seq_df, file = "otu_sequences_97percent.csv", row.names = FALSE)



### 6. OTU curation ####

# File preparation
# A - OTU table with samples as columns and OTUs as rows
tax_dada_otu  = t(tax_dada_otu)
tax_dada_otu.df  = phyloseq_to_df(tax_dada_otu, addtax = F, addtot = F, addmaxrank = F,
                                  sorting = "abundance") # with package metagMisc
tax_dada_otu.df  <- data.frame(tax_dada_otu.df , row.names = 1)

# B - Fasta file to prepare the match list
# from manually created .csv file
fa.table = read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\otu_sequences_97percent.csv")
fa = character(2 * nrow(fa.table))
fa[c(TRUE, FALSE)] = paste0(">", fa.table$ASV)
fa[c(FALSE, TRUE)] = as.character(fa.table$Sequence)

writeLines(fa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\otus97.fasta")

# C - Create match list with vsearch
# Place input file "otus97.fasta" in "/Applications/vsearch/bin" folder of vsearch and run:
# vsearch --usearch_global otus97.fasta --db otus97.fasta --self --id .84 --iddef 1 --userout match_list97.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

matchlist_name = read.table("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\vsearch-2.30.0-win-x86_64\\vsearch-2.30.0-win-x86_64\\bin\\match_list.txt")
names(matchlist_name)[names(matchlist_name) == "V1"] <- "OTUid"
names(matchlist_name)[names(matchlist_name) == "V2"] <- "hit"
names(matchlist_name)[names(matchlist_name) == "V3"] <- "match"
matchlist_name$OTUid <- as.character(matchlist_name$OTUid)
matchlist_name$hit <- as.character(matchlist_name$hit)


# Run OTU curation
curated_result <- lulu(tax_dada_otu.df, matchlist_name)

# Curated OTU table
curated_table = curated_result$curated_table
write.csv(curated_table, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\curated_table97.csv")

# Curation Results 
# 550 97% OTUS were obtained. 

# Prepare fasta file for taxonomic assignment using the curated table. 
# Create sequence list by saving first column of "curated_table97.csv" as a new file.
# I called mine "curated97_seq.csv". 
##Como este meu arquivo contém apenas os nomes das sequencias, adicionei uma parte do código que associasse os nomes que sobraram no ultimo passo com as sequencias as quais eles correspondem

# Ler os dois arquivos
asv_names <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\curated97_names.csv", stringsAsFactors = FALSE)
asv_seqs <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\Alcatrazes_2023_raw\\filtered\\Dereplicated\\taxonomic_assignment\\5_Clustering_ASV_into_OTU\\otu_sequences_97percent.csv", stringsAsFactors = FALSE)

# Verifique os nomes das colunas para garantir que ambas tenham "asv"
# names(asv_names)
# names(asv_seqs)

# Fazer o merge usando a coluna 'asv'
merged_df <- merge(asv_names, asv_seqs, by = "ASV")

# Salvar o novo arquivo
write.csv(merged_df, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\asv_com_sequencias.csv", row.names = FALSE)

##

curated_seqs = read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\asv_com_sequencias.csv")
fa = character(2 * nrow(curated_seqs))
fa[c(TRUE, FALSE)] = paste0(">", curated_seqs$ASV)
fa[c(FALSE, TRUE)] = as.character(curated_seqs$Sequence)
writeLines(fa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\curated_97OTUseqs_its2.fasta") # name typo


### 7. Repeat taxonomic assignment with insect ####
# --> Done by Helio Quintero at NAOS Computers. 
# --> The classifier used was provided by insect github
# --> https://github.com/shaunpwilkinson/insect

# Script used:

x <- readFASTA("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\curated_97OTUseqs_its2.fasta")
classifier <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\classifier.rds") 
out <- classify(x, classifier, ping = 0.98, cores = 8)

write.csv(out, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\otus97_insect.csv")

### Note: I rename this csv as "otus97_insect.csv" in my laptop. 

#_________________________________________________________________________
#
#### 8. New Phyloseq ####
#_________________________________________________________________________


# Transpose otu table
otu_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\curated_table97.csv")
asv_mat <- t(otu_mat)
write.csv(asv_mat, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\otumat97.csv")

# Read tables
otu_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\otumat97_v2.csv")
tax_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\otus97_insect_v2.csv")
sample_df <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\SAM_otu.csv")

# Define rownames
otu_mat <- otu_mat %>% 
  tibble::column_to_rownames("X") #era "asv"
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("representative") #era "asv"
sample_df <- sample_df %>%
  tibble::column_to_rownames("SampleID") #era "sample"

# Transpõe para que ASVs fiquem nas linhas
otu_mat <- t(otu_mat)

# Convert to matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat) 

# Create tabs 
asv = otu_table(otu_mat, taxa_are_rows = T)
tax = tax_table(tax_mat)
sample = sample_data(sample_df)

# Create phyloseq
alcatrazes_2023_otu <- phyloseq(asv, tax, sample)
alcatrazes_2023_otu

head(taxa_names(asv))
head(taxa_names(tax))

saveRDS(alcatrazes_2023_otu, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\alcatrazes_2023_otu.rds")

#Verificar estrutura geral

summary(alcatrazes_2023_otu)


#Verificar riqueza/diversidade

plot_richness(alcatrazes_2023_otu, x = "Site", measures = c("Shannon", "Observed")) 

#Checar composição por filo, etc...

unique(tax_table(alcatrazes_2023_otu)[,"phylum"])
relative_filtered <- subset_taxa(alcatrazes_2023_otu, phylum != "")
plot_bar(relative_filtered, fill = "phylum")

# Exclui OTUs dos filos Ctenophora e Nematoda
alcatrazes_filtrado <- subset_taxa(
  alcatrazes_2023_otu,
  !phylum %in% c("Ctenophora", "Nematoda")
)


#Verificar riqueza/diversidade

plot_richness(alcatrazes_filtrado, x = "Site", measures = c("Shannon", "Observed")) #Substitua "SEU_GRUPO" pelo nome de uma variável categórica em sample_data, como Site, Depth, Treatment, etc.


#Checar composição por filo, etc...

unique(tax_table(alcatrazes_filtrado)[,"phylum"])
relative_filtered <- subset_taxa(alcatrazes_filtrado, phylum != "")
plot_bar(relative_filtered, fill = "phylum")


#### CR Report 

# Get on the working directory 
setwd("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq")

rohr42 <- readRDS("alcatrazes_2023_otu.rds") #carregou objeto phyloseq
rohr42  <- prune_taxa(taxa_sums(rohr42) > 2, rohr42) #removeu OTUs com menos de 3 leituras
rohr42.nona = subset_taxa(rohr42, !is.na(phylum) & !phylum %in% c("NA", "")) #removeu taxons sem identificação de filo
rohr.norm = transform_sample_counts(rohr42.nona, function(x) x / sum(x) * 100) #transformou abundâncias em porcentagens
rohr42.target <- subset_taxa(rohr42.nona, phylum %in% c("Porifera", "Cnidaria")) #manter apenas dados de poriferos e cnidarios
rohr42.target.norm <- transform_sample_counts(rohr42.target, function(x) x / sum(x) * 100) #transformar em porcentagem


#Esta etapa é para remover outliers, n fiz pq n sei oq eu poderia considerar outlier 
# Open new phyloseq with the classification of plancto and bentho
#setwd("C:/Users/Viviane/Documents/2024/Benthic_eDNA")
#rohr <- readRDS("coralotus2.rds")
# 8.1 Take out outliers   ####
#rohr_out <- subset_samples(rohr,
#                           sample_names(rohr) != "ML3997_ITS2_SABA"&
#                             sample_names(rohr) != "ML3635_ITS2_BAD" &
#                             sample_names(rohr) != "ML3636_ITS2_BAD" &
#                             sample_names(rohr) != "ML3910_ITS2_WS" &
#                             sample_names(rohr) != "ML3887_ITS2_PAJ" & 
#                             sample_names(rohr) != "ML3909_ITS2_WS" &
#                             sample_names(rohr) != "ML3911_ITS2_WS" & 
#                             sample_names(rohr) != "ML3638_ITS2_BAD" &
#                             sample_names(rohr) != "ML3889_ITS2_PAJ" & 
#                             sample_names(rohr) != "ML3888_ITS2_PAJ")
#saveRDS(rohr_out, "C:/Users/Viviane/Documents/2024/Benthic_eDNA/rohr_out.rds")


# 8.2 Relative abundance    ####

relative <- rohr42.target.norm

#Aplicar filtros
# Remove amostras com valores ausentes
relative <- prune_samples(!is.na(sample_sums(relative)), relative)
# Remove táxons com valores ausentes
relative <- prune_taxa(!is.na(taxa_sums(relative)), relative)
#verificando se funcionou, deve retornar FALSE
any(is.na(otu_table(relative)))

saveRDS(relative, file = "phyloseq_filtrado_relative.rds") #este é um novo arquivo phyloseq após o inicial ser filtrado, como eu ja chamei ele de "relative" nao precisei colocar o caminho para ele nas proximas etapas. mas ainda assim estamos lidando com os menos dados pq "relative" faz referencia ao phyloseq inicial filtrado :)
#antes era rohr_out no script da viviane


df <- sample_data(relative) %>% as.data.frame()

variable.names(df)

nmds <- ordinate(relative, "NMDS", "bray")
nmds_plot <- plot_ordination(relative, nmds, color = "Site", shape = "Season")
nmds_plot + 
  geom_point(size = 4) + 
  scale_color_manual(values = c("tomato3", "turquoise4", "paleturquoise3", "purple", "orange", "pink"),
                     labels = list(
                       "P1",
                       "P2",
                       "P3",
                       "P4",
                       "P5",
                       "P6"
                     )) + 
  ggtitle("Composição porifera e cnidaria alcatrazes por estação em 2023") +
  theme_minimal() + 
  stat_ellipse(linewidth = 1.5, alpha = 0.3, level = 0.98) + 
  theme(legend.title = element_blank()) + 
  annotate(geom = "text", 
           label = "PERMANOVA\nF = 26.45, p-value = 0.001", 
           x = -1.5, y = -1, size = 3, family = "sans")
nmds_plot

#A seguinte etapa deu erro e eu modifiquei
#df <- data.frame(sample_data(relative))
#dist <- distance(relative, method = "bray")
#perma <- adonis2(dist ~ Site, df)
#stress <- metaMDS(dist)
#disper <- betadisper(dist, df$Site, type = "centroid")
#permt <- permutest(disper, pairwise = T)
#permt


#Fiz assim:
# Extrair a tabela de OTUs do objeto phyloseq
otu_table_relative <- otu_table(relative)

# Verificar se a tabela OTU tem a estrutura correta
head(otu_table_relative)

# Calcular a distância de Bray-Curtis
dist <- phyloseq::distance(relative, method = "bray")

# Análise de PERMANOVA
df <- data.frame(sample_data(relative))
perma <- adonis2(dist ~ Site, df)
print(perma)

# Calcular o stress para NMDS
stress <- metaMDS(dist)
print(stress)

# Calcular a dispersão beta (Betadisper)
disper <- betadisper(dist, df$Site, type = "centroid")
print(disper)

# Rodar o teste de permutação para dispersão
permt <- permutest(disper, pairwise = T)
print(permt)


# by Site 
# Se você quer usar o objeto 'relative' que já está filtrado e normalizado:
CR <- relative

# NMDS por Site
nmds <- ordinate(CR, method = "NMDS", distance = "bray")

nmds_plot <- plot_ordination(CR, nmds, color = "Site") +
  geom_point(size = 4) +
  scale_color_manual(values = c("tomato3", "tomato4", "navajowhite4", "olivedrab",
                                "cornsilk3", "burlywood3")) +
  ggtitle("Poriferos e cnidarios alcatrazes") +
  theme_minimal() +
  theme(legend.title = element_blank())

nmds_plot

#permanova
dfcr <- data.frame(sample_data(CR))

# Usar distância Bray-Curtis com vegan::vegdist
library(vegan)
otu <- as.data.frame(otu_table(CR))
if (taxa_are_rows(CR)) otu <- t(otu)

distcr <- vegdist(otu, method = "bray")

# PERMANOVA
permacr <- adonis2(distcr ~ Site, data = dfcr)
print(permacr)

# Análise par a par (certifique-se de ter o pacote disponível)
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)

pairwise.adonis2(distcr ~ Site, data = dfcr)


#Plotar a dispersão multivariada dos grupos (Site):
# Análise de dispersão
disper <- betadisper(distcr, dfcr$Site)
# Teste de permutação para homogeneidade dos grupos
permutest(disper)
# Plot com centroides e elipses
plot(disper, hull = FALSE, ellipse = TRUE, label = TRUE, main = "Dispersão por Site")


# 8.3 Rarefaction curves to see the depth of sequencing    ####

devtools::install_github("gauravsk/ranacapa")
library(ranacapa)

sample_sums(rohr42)
rohr42 <- prune_samples(sample_sums(rohr42) >= 500, rohr42)

p0 <- ggrare(rohr42,
             step = 500,
             color = "Site",
             plot = T,
             parallel = T,
             se = F)
p0 <- p0 + 
  facet_wrap(~ Site) + 
  geom_line(linewidth=0.8) +
  geom_vline(xintercept = min(sample_sums(rohr42)), 
             color = "gray60") +
  theme_bw()

plot(p0)

rohr.pa <- transform_sample_counts(rohr42, function(abund)1*(abund>0))
richness <- as.data.frame(sample_sums(rohr.pa))

ggplot(richness, aes(`sample_sums(rohr.pa)`)) +
  geom_histogram() +
  labs(
    title = "Samples Histogram", 
    subtitle = "Observed species in each sample: cnidaria and porifera eDNA", 
    caption = "data from alcatrazes 2023"
  )


# Take out the outliers 

#planctonic <- subset_taxa(rohr_out, type == "planthonic") #"rohr_out" seria o "relative", mas nao temos essa diferenciação de planc e bent no nosso phyloseq
#benthonic <- subset_taxa(rohr_out, type == "benthonic")

sclerac <- subset_taxa(relative, order == "Scleractinia")
taxa.scler <- as.data.frame(sclerac@tax_table)
write.csv(taxa.scler, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\Alcatrazes_2023\\6_OTU_curation\\7_repeat_taxonomic_assignment\\8_new_phyloseq\\scleractinia_taxtable.csv")

#de novo, n fizemos essa separação, ent n fiz isso
#bentho_rabund <- transform_sample_counts(benthonic, function(x) x/sum(x) *100) 
#plantho_rabund <- transform_sample_counts(planctonic, function(x) x/sum(x) *100)


#########More exploratory graphs##############



# Barras empilhadas por amostra por ordem, etc...
unique(tax_table(relative)[,"order"])
relative_filtered <- subset_taxa(relative, order != "")
plot_bar(relative_filtered, fill = "order") 

# Barras empilhadas com cores definidas por ordem, etc...
plot_bar(relative, fill = "order") + 
  geom_bar(aes(color = order, fill = order), stat = "identity", position = "stack")

# Facetando por variável do metadado para ordem, etc...

unique(tax_table(relative)[,"order"])
relative_filtered <- subset_taxa(relative, order != "")
plot_bar(relative_filtered, x = "order", facet_grid = ~Site, fill = "order") +
  scale_fill_manual(values = c(
    "Zoantharia" = "blue",
    "Verongiida" = "red",
    "Scleractinia" = "green",
    "Leptothecata" = "pink",
    "Ceriantharia" = "orange",
    "Alcyonacea" = "purple",
    "Actiniaria" = "yellow"
  ))


unique(tax_table(relative)[,"family"])
relative_filtered <- subset_taxa(relative, family != "")
plot_bar(relative_filtered, x = "family", facet_grid = ~Site, fill = "family") 


unique(tax_table(relative)[,"genus"])
relative_filtered <- subset_taxa(relative, genus != "")
plot_bar(relative_filtered, x = "genus", facet_grid = ~Site, fill = "genus") 


unique(tax_table(relative)[,"genus"])
relative_filtered <- subset_taxa(relative, genus != "")
plot_bar(relative_filtered, x = "family", facet_grid = ~Season, fill = "genus")



############################################

# 8.4 betapart analysis  ####

library("betapart")

sitio <- c("Coiba", "Perlas", "Cocos")

beta_tab <- read.csv("betapart_ITS2.csv", 
                     row.names = sitio)

betapart <- decostand(beta_tab[,-1], 
                      method = "pa")

core <- betapart.core(betapart)
incidence <- beta.multi(core, index.family = "sorensen")
beta_samp <- beta.sample(core, 
                         index.family = "sor", 
                         sites = 3, 
                         samples = 1000)


# Plotting 

plot(density(beta_samp$sampled.values$beta.SIM), col = "red", 
     xlim = c(0,0.9), ylim = c(0,14), lwd = 1.5)
lines(density(beta_samp$sampled.values$beta.SNE), col = "blue", 
      xlim = c(0, 0.9), ylim = c(0,14), lwd = 1.5)
lines(density(beta_samp$sampled.values$beta.SOR), col = "black", 
      xlim = c(0, 0.9), ylim = c(0,14), lwd = 1.5)


# Extracting the otu tables of benthic and plankto 

otu_benth <- as.data.frame(benthonic@otu_table)
write.csv(otu_benth, "C:/Users/Viviane/Documents/2024/Benthic_eDNA/otu_benth.csv")

otu_plant <- as.data.frame(planctonic@otu_table)
write.csv(otu_plant, "C:/Users/Viviane/Documents/2024/Benthic_eDNA/otu_plant.csv")



# 8.5 NMDS by sites ####

perlas <- subset_samples(rohr_out, rohr_out@sam_data$Locality == "Las Perlas")
coiba <- subset_samples(rohr_out, rohr_out@sam_data$Locality == "Coiba island")
cr <- subset_samples(rohr_out, rohr_out@sam_data$Locality == "Cocos island")

perlas.rabund <- transform_sample_counts(perlas, function(x) x/sum(x)*100)
coiba.rabund <- transform_sample_counts(coiba, function(x) x/sum(x)*100)
cr.rabund <- transform_sample_counts(cr, function(x) x/sum(x)*100)


nmds_perlas <- ordinate(perlas.rabund, "NMDS", "bray") %>% 
  plot_ordination(perlas.rabund, ., color = "site")
nmds_perlas <- nmds_perlas + 
  geom_point(size = 4) +
  scale_color_manual(values = c("paleturquoise3", "lightcyan3", "honeydew2")) + 
  ggtitle("Composición Bentónica") +
  stat_ellipse(linewidth = 1.5, alpha = 0.3, level = 0.98) +
  theme_minimal() + 
  theme(legend.title = element_blank())

nmds_perlas

df_perlas <- data.frame(sample_data(perlas.rabund))
dist_perlas <- distance(perlas.rabund, method = "bray")
perma_perlas <- adonis2(dist_perlas ~ site, df_perlas)
pair_perlas <- pairwise.adonis2(dist_perlas ~ site, df_perlas)



nmds_coiba <- ordinate(coiba.rabund, "NMDS", "bray") %>% 
  plot_ordination(coiba.rabund, ., color = "site")
nmds_coiba <- nmds_coiba + 
  geom_point(size = 4) +
  scale_color_manual(values = c("turquoise4", "lightseagreen", "mediumaquamarine")) + 
  ggtitle("Composición Bentónica") +
  theme_minimal() + 
  theme(legend.title = element_blank())

nmds_coiba


df_coiba <- data.frame(sample_data(coiba.rabund))
dist_coiba <- distance(coiba.rabund, method = "bray")
perma_coiba <- adonis2(dist_coiba ~ site, df_coiba)
pair_coiba <- pairwise.adonis2(dist_coiba ~ site, df_coiba)



nmds_cr <- ordinate(cr.rabund, "NMDS", "bray") %>% 
  plot_ordination(cr.rabund, ., color = "site")
nmds_cr <- nmds_cr + 
  geom_point(size = 4) +
  scale_color_manual(values = c("tomato3","tomato4", 
                                "navajowhite4","cornsilk3")) + 
  ggtitle("Composición Bentónica") +
  stat_ellipse(linewidth = 1.5, alpha = 0.3, level = 0.98) +
  theme_minimal() + 
  theme(legend.title = element_blank())

nmds_cr


df_cr <- data.frame(sample_data(cr.rabund))
dist_cr <- distance(cr.rabund, method = "bray")
perma_cr <- adonis2(dist_cr ~ site, df_cr)
pair_cr <- pairwise.adonis2(dist_cr ~ site, df_cr)



nmds_bentho <- ordinate(bentho_rabund, "NMDS", "bray") %>% 
  plot_ordination(bentho_rabund, ., color = "Locality", shape = "season")
nmds_bentho <- nmds_bentho + 
  geom_point(size = 4) +
  scale_color_manual(values = c("tomato3","turquoise4", "paleturquoise3"),
                     labels = list(
                       "Isla del Coco, CR",
                       "Isla Coiba, PA",
                       "Archipiélago de Las Perlas, PA"
                     )) + 
  ggtitle("Composición Bentónica") +
  theme_minimal() + 
  stat_ellipse(linewidth = 1.5, alpha = 0.3, level = 0.98) + 
  theme(legend.title = element_blank())
nmds_bentho



nmds_plancto <- ordinate(plantho_rabund, "NMDS", "bray") %>% 
  plot_ordination(plantho_rabund, ., color = "Locality", shape = "season")
nmds_plancto <- nmds_plancto + 
  geom_point(size = 4) +
  scale_color_manual(values = c("tomato3","turquoise4", "paleturquoise3"),
                     labels = list(
                       "Isla del Coco, CR",
                       "Isla Coiba, PA",
                       "Archipiélago de Las Perlas, PA"
                     )) + 
  ggtitle("Composición Plantónica") +
  theme_minimal() + 
  stat_ellipse(linewidth = 1.5, alpha = 0.3, level = 0.98) + 
  theme(legend.title = element_blank())
nmds_plancto



# 8.6 Extract cnidarian ASVs ####

rohr

cnidarians <- subset_taxa(rohr, phylum == "Cnidaria")
cnida <- tax_glom(cnidarians, taxrank = "taxon")

cnida_df <- as.data.frame(cnida@tax_table)

write.csv(cnida_df, "C:/Users/Viviane/Documents/2024/Benthic_eDNA/cnida_df.csv")


