# https://alexiscarter.github.io/metab/

# https://alexiscarter.github.io/metab/Dada_script_EN.html

# https://benjjneb.github.io/dada2/tutorial.html


#12S

#https://lizsuter.github.io/files/DADA2_pipeline_SCM_eDNA.nb.html


######################################################### Install and load packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("phyloseq")

install.packages(c('ggplot2', 'vegan', 'gtools'))
install.packages('reshape2')

library(dada2); packageVersion("dada2")
library (vegan)
library (gtools)
library (ggplot2)
library (phyloseq)
library(reshape2)

########################################################### set work directory

setwd("~/Desktop/fish_12s_2023/raw_data/12S_Lib2/12S_Lib2_demultiplexado/")

########################################################### check list files, set path, and rename samples

list.files()

path <- "~/Desktop/fish_12s_2023/raw_data/12S_Lib2/"
fnFs <- sort(list.files(pattern="_R1.fastq.gz"))
fnRs <- sort(list.files(pattern="_R2.fastq.gz"))

sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names

#plotQualityProfile(fnFs[1:4]) # use essa funcao caso queira restringir a analise
plotQualityProfile(fnFs)
############################################################ Quality profile

fnFs_plotQP <- list()
fnRs_plotQP <- list()

for (i in 1:74) {
  pltName <- sample.names[i]
  fnFs_plotQP[[pltName]] <- plotQualityProfile(fnFs[(i)])
  fnRs_plotQP[[pltName]] <- plotQualityProfile(fnRs[(i)])
}

#### Arrange the plots in a single image (With 3 samples)

library(ggpubr)

fnFs_arrange <- ggarrange(fnFs_plotQP$PC_fracao1_Tag2_R1.fastq, fnFs_plotQP$PC_fracao2_Tag3_R1.fastq, fnFs_plotQP$PC_fracao3_Tag5_R1.fastq, 
                          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                          align = "hv",
                          widths = 1,
                          heights = 1,
                          label.x = 0,
                          label.y = 1,
                          hjust = -0.5,
                          vjust = 1.5,
                          ncol = 2, nrow = 2)

fnFs_arrange

fnRs_arrange <- ggarrange(fnRs_plotQP$PC_fracao1_Tag2_R1.fastq, fnRs_plotQP$PC_fracao2_Tag3_R1.fastq, fnRs_plotQP$PC_fracao3_Tag5_R1.fastq, 
                          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                          align = "hv",
                          widths = 1,
                          heights = 1,
                          label.x = 0,
                          label.y = 1,
                          hjust = -0.5,
                          vjust = 1.5,
                          ncol = 2, nrow = 2)
fnRs_arrange

fnFs_arrange
ggsave("fnFs_arrange.png", device = "png", width = 30, height = 20, units = "cm", dpi = 300)

fnRs_arrange
ggsave("fnRs_arrange.png", device = "png", width = 30, height = 20, units = "cm", dpi = 600)

############################################################ Filtering and Trimming

# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered_pairedend") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#run filter and trim sequences - 515FY 19BASES, 929 20BASES

#primerForward_Tarbelet_et_al_2018_Teleo02 <-"XXXXXXXXAAACTCGTGCCAGCCACC"
#primerReverse_Tarbelet_et_al_2018_Teleo02 <- "XXXXXXXXGGGTATCTAATCCCAGTTTG"

#explaining parameters
#truncQ=2 #PHRED quality score # https://forum.qiime2.org/t/why-do-the-dada2-default-setting-have-such-a-low-phred-score/10156
#trimLeft=c(32,32) #used to specify the number of bases to trim from the beginning (5' end) of each read before quality filtering and truncation # 32 means lenght of primers
#maxEE=c(2,2)) #stands for the maximum expected error rate, which is used as a quality control threshold during the denoising step # The value of 2 is often used as a conservative threshold, meaning that any sequences with an expected error rate greater than 2 will be discarded.
#trunclen # length for the forward and reverse reads after filter and trimming # By adjusting the truncLen parameter, you can control the length of the sequences used for downstream analysis in DADA2
#just in case you need to remove primers/adapter from cutadapt # https://benjjneb.github.io/dada2/ITS_workflow.html

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncQ=2, #default dada2
                     trimLeft=c(26,24), #tamanho dos primers (26pb) + tag (6pb) (mlCOIintF: XXXXXXGGWACWGGWTGAACWGTWTAYCCYCC e jgHCO2198: XXXXXXTAIACYTCIGGRTGICCRAARAAYCA - X means 6 tags)
                     maxEE=c(2,2))  #default dada2
#multithread=TRUE)
head(out)

############################################################ Filtering visualization

pourc <- cbind((out[,2]/out[,1])*100) # Percentage filtered sequence / non-filtered sequence
pourc_disc <- cbind(out, pourc) # combines out and pourc
pourc_disc 

(mean(out[,2])/mean(out[,1]))*100 # Mean percentage 

############# Create list of quality plots and iterate over the samples
fnFs_filt_plotQP <- list()
fnRs_filt_plotQP <- list()

#plot das sequências filtradas e trimadas
plotQualityProfile(filtFs[1]) # 1st Forward sample
plotQualityProfile(filtRs[1]) # 1st Reverse sample
plotQualityProfile(filtFs, aggregate = TRUE)
plotQualityProfile(filtRs, aggregate = TRUE)


for (i in 1:4) {
  pltName <- sample.names[i]
  fnFs_filt_plotQP[[pltName]] <- plotQualityProfile(filtFs[(i)])
  fnRs_filt_plotQP[[pltName]] <- plotQualityProfile(filtRs[(i)])
}

##### Arrange the plots in a single image #####


fnFs_filt_arrange <- ggarrange(fnFs_filt_plotQP$PC_fracao1_Tag2_R1.fastq, fnFs_filt_plotQP$PC_fracao2_Tag3_R1.fastq, fnFs_filt_plotQP$PC_fracao3_Tag5_R1.fastq,
                               font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                               align = "hv",
                               widths = 1,
                               heights = 1,
                               label.x = 0,
                               label.y = 1,
                               hjust = -0.5,
                               vjust = 1.5,
                               ncol = 2, nrow = 2)

fnFs_filt_arrange

ggsave("fnFs_filt_arrange.png", device = "png", width = 30, height = 20, units = "cm", dpi = 300)


fnRs_filt_arrange <- ggarrange(fnRs_filt_plotQP$PC_fracao1_Tag2_R1.fastq, fnRs_filt_plotQP$PC_fracao2_Tag3_R1.fastq, fnRs_filt_plotQP$PC_fracao3_Tag5_R1.fastq, 
                               font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                               align = "hv",
                               widths = 1,
                               heights = 1,
                               label.x = 0,
                               label.y = 1,
                               hjust = -0.5,
                               vjust = 1.5,
                               ncol = 2, nrow = 2)

fnRs_filt_arrange


ggsave("fnRs_filt_arrange.png", device = "png", width = 30, height = 20, units = "cm", dpi = 300)


#sites ro help interprated these type of plot
#https://github.com/benjjneb/dada2/issues/1307
#https://github.com/benjjneb/dada2/issues/1117


############################################################# Error Rates Learning

# This step consist in estimating the error rates due to sequencing. Its purpose is to differentiate between mutant sequences and false sequences.

errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

tiff("errF.tif", width = 10, height = 10, units = "in", res = 300)
plotErrors(errF, nominalQ=TRUE)
dev.off()

tiff("errR.tif", width = 10, height = 10, units = "in", res = 300)
plotErrors(errR, nominalQ=TRUE)
dev.off()

save(errF, errR, file = "dada2_after_errF_errR.RData")

############################################################ Dereplicating
# Combines all identical sequencing reads into into unique sequences with a corresponding abundance
# perform duplicate removal on a filtered Fastq file
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names

derepRs <- derepFastq(filtRs)
names(derepRs) <- sample.names

############################################################ Sample Inference
dadaFs <- dada(derepFs, 
               err = errF, 
               multithread=TRUE,
               pool=TRUE)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, 
               err=errR,
               multithread=TRUE,
               pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

save(dadaRs, file="dadaRs.rdata")
save(dadaFs, file="dadaFs.rdata")

########################################################### Merging paired reads

#The minOverlap parameter is typically set to an integer value that represents the minimum number of overlapping bases required for merging. 
#By default, DADA2 sets minOverlap to 12.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 12, 
                      maxMismatch = 0)

head(mergers[[1]])
max(mergers[[1]]$nmatch) # Largest overlap # it does not work
min(mergers[[1]]$nmatch) # Smallest overlap #it does not work     

########################################################### ASVs table
#construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab[,1]
seqtab

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)),xlab="Size", ylab="Frequency", main = "ASVs length", xlim=c(100,300), ylim=c(0,5000)) 

########################################################### Removing chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "per-sample", 
                                    multithread = TRUE,
                                    verbose = TRUE)

sum(seqtab.nochim)/sum(seqtab)


save(seqtab.nochim, file="seqtab.nochim.rdata")


#load("seqtab.nochim.rdata")

########################################################### inspect the results

round((sum(seqtab.nochim)/sum(seqtab)*100),2) # Percentage of the total sequence reads

hist(nchar(getSequences(seqtab.nochim)),xlab="Size", ylab="Frequency", main = "Non-chimeric ASVs length", xlim=c(100,250), ylim=c(0,2500)) # Lenght of the non-chimeric sequences

########################################################### transform the ASVs occurrences in presence / absence 

# presence / absence allows to quantify the number of ASVs per sample.

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) 

########################################################## track Table - it doesn't work when using just one sample

getN <- function(x) sum(getUniques(x))
track1 <- cbind(out, sapply(dadaFs,getN),sapply(dadaRs,getN),sapply(mergers,getN),rowSums(seqtab.nochim))
colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track1) <- sample.names
head(track1)

write.csv(track1, "./track_table.csv")

## figure overall analysis
#do it just if you have more than one sample

## figure overall analysis
#do it just if you have more than one sample

gtrack1<- track1[,c(2-7)]
gtrack1$ID <- rownames()

lgtrack1 <- melt(gtrack1, id.vars="ID")
View(lgtrack1)
bar_track1 <- ggplot(lgtrack1 ,aes(x=Var1, y=as.numeric(value), fill=Var2)) +
  geom_bar(stat="identity", position = "identity") + 
  theme_classic() + # Theme
  theme(axis.ticks.length=unit(0.3,"cm")) + # Ticks size
  theme(axis.text.x = element_text(angle=45) , legend.title = element_blank())+ # Changes the x labels orientation & delete legend title
  scale_x_discrete(name ="Sample ID", limits=rownames(track1))+ # Changes x-axis title & sorts the x label names
  scale_y_continuous(name="Abundance", breaks=seq(from = 0, to = 1000, by = 100))+ #Changes y-axis title & sets the y breaks.
  ggtitle("Track")# Main title
bar_track1

ggsave("seq_track1.png", device = "png", width = 30, height = 20, units = "cm", dpi = 300)

head (seqtab.nochim)

###### Links to continue processing 
### https://github.com/elaine-shen/Indo_eDNA_RSMS/blob/main/4.%20Two-part%20taxonomic%20assignment%20using%201)%20CO1v4%20database%20and%20RDP%20Classifier%20and%202)%20BASTA/combine_taxassign.R
### https://github.com/terrimporter/CO1Classifier?tab=readme-ov-file

## Create a fasta file for the unique sequences and save final combined sequence table

uniquesToFasta(getUniques(seqtab.nochim), fout="~/Desktop/fish_12s_2023/raw_data/12S_Lib2/uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))
saveRDS(seqtab.nochim, "~/Desktop/fish_12s_2023/raw_data/12S_Lib2/combined.rds") # CHANGE ME to where you want sequence table saved

#### MAKE ASV FASTA FILE FOR TAXONOMIC ASSIGNMENT ####
#seqtab_nochim = readRDS("/home/lagem/Documentos/ARMS_Lucas/2024/CC_S738/combined.rds")

sq <- getSequences(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}
names(sq) <- asv_headers

library(ShortRead)
writeFasta(sq, file="combined-asvs.fasta")

# Make taxonomy table with ASVID and sequence 
library (devtools)
library (tidyverse)

source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
FastaToTabular("combined-asvs.fasta")

############################################################################# CONTINUAR O SCRIPT


#################################
##### CURATE ASVS WITH LULU #####
#################################

# https://github.com/elaine-shen/Indo_eDNA_RSMS/blob/main/3.%20Create%20OTU%20table%20from%20ASVs/LULU_ASVs.R

library(devtools)
install_github("tobiasgf/lulu")  
library(lulu)

# PREP FOR LULU 
library(openssl)


# Functions ------------------------------------
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(getwd(), "DADA2_extracted_samples")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",openssl::sha1(colnames(my_table)[seqX]),";size=",
                           my_table[sampleX,seqX],";barcodelabel=",rownames(my_table)[sampleX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}


# generate tab-delimited OTU table with taxa as rows
DADA2_otutab <- readRDS("/home/lagem/Documentos/ARMS_Lucas/2024/PC_S738/combined.rds")
LULU_DADA2_otutab <- t(DADA2_otutab)
# make rownames SHA1 sums
print(rownames(LULU_DADA2_otutab)[1:5])
rownames(LULU_DADA2_otutab) <- sha1(rownames(LULU_DADA2_otutab))
write.table(LULU_DADA2_otutab, file = "DADA2_sha1_otutab.tsv", sep = "\t", quote = FALSE)

extrSamDADA2(DADA2_otutab)


#############################
##### combine taxassign #####
#############################


setwd("/Users/elaineshen/Desktop/cox1_raw/combined/databases/CO1v4_training/mydata_trained")

##### LIBRARIES ######
library(dplyr)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(microbiome) # accessing more of phyloseq
library(vegan) # rarefaction curves
library(metagMisc) # format for iNEXT
library(iNEXT) # rarefaction curves


##### IMPORT RAW DATA ######
## using excel sheets here instead of text, has the same effect ##

## Taxonomy files
setwd("~E:/Metabarcoding_ARMS_SP/2024/CA_S738")
seqs_raw = readxl::read_xlsx("E:/Metabarcoding_ARMS_SP/2024/CA_S738//COIv5.output.xlsx")
colnames(seqs_raw) = c("name","V2", "root", "root2","root_confidence","domain_name","domain_rank", "domain_confidence", 
                       "kingdom_name", "kingdom_rank", "kingdom_confidence",
                       "phylum_name", "phylum_rank", "phylum_confidence",
                       "class_name", "class_rank","class_confidence",
                       "order_name", "order_rank", "order_confidence",
                       "family_name", "family_rank","family_confidence",
                       "genus_name", "genus_rank", "genus_confidence",
                       "species_name", "species_rank","species_confidence")

## avaluate necessity of BASTA search
# Iterative BLAST output with LCA 
ib_output = read.delim("/Users/elaineshen/Desktop/cox1_raw/combined/iterative_blast/ib_output/all_blast_basta.txt", header=F)
ib_output = ib_output[!duplicated(ib_output[c(1,2)]),] # Remove duplicate assignments for an ASV if they are the same
ib_output = ib_output %>% filter(!str_detect(V2, 'Unknown')) # Get rid of "Unknown" IB outputs
colnames(ib_output) = c("name", "tax")

ib_output = separate(data = ib_output, col = "tax", into = c("domain_name", "phylum_name","class_name","order_name","family_name","genus_name","species_name"), sep = "\\;")
ib_output[ ,c("domain_confidence","phylum_confidence","class_confidence","order_confidence", "family_confidence","genus_confidence","species_confidence")] <- 0.9701010101 # Minimum possible confidence for the iterative blast and the LCA was 97, so set all of these rank confidences as 0.9701010101 to separate these from RDP outputs
ib_output$domain_rank = "superkingdom"
ib_output$phylum_rank = "phylum"
ib_output$class_rank = "class"
ib_output$order_rank = "order"
ib_output$family_rank = "family"
ib_output$genus_rank = "genus"
ib_output$species_rank = "species"

# Add kingdom info as undef_domain
ib_output$kingdom_name = paste0("undef_",ib_output$domain_name)
ib_output$kingdom_rank = "kingdom"
ib_output$kingdom_confidence = 0.9701010101 # This output did not provide kingdom-level information, only domain. But will assume 97% to be conservative

# Add misc columns for ease of merging
ib_output[,c("root","root2")] = "cellularOrganisms"
ib_output$V2 = ""
ib_output$root_confidence = 1 # all should be 1

# Reorder columns
col_order = c("name","V2", "root", "root2","root_confidence","domain_name","domain_rank", "domain_confidence", 
              "kingdom_name", "kingdom_rank", "kingdom_confidence",
              "phylum_name", "phylum_rank", "phylum_confidence",
              "class_name", "class_rank","class_confidence",
              "order_name", "order_rank", "order_confidence",
              "family_name", "family_rank","family_confidence",
              "genus_name", "genus_rank", "genus_confidence",
              "species_name", "species_rank","species_confidence")
ib_output = ib_output[ ,col_order]


# Replace the CO1v4 taxonomy outputs with the iterative blast output
seqs_raw[match(ib_output$name,seqs_raw$name), ] <- ib_output

## ASV table
asv_table = as.data.frame(readRDS("/Users/elaineshen/Desktop/cox1_raw/combined/combined.rds"))

asv_tax = read.csv("/Users/elaineshen/Desktop/cox1_raw/combined/asv_table.csv", header=TRUE)
asv_tax$name<-gsub(">","",as.character(asv_tax$name))
asv_tax = asv_tax[,c(2,3)]

# Need to replace asv table's columns, which are the sequences, with the ASV ID. I could also add the sequence to the taxonomy table (but that's less clean)
names(asv_table)[match(asv_tax[,"sequence"], names(asv_table))] = asv_tax[,"name"]



##### CREATING GLOBAL CONFIDENCE/TAXONOMY COLUMN, REMOVING IDs WITH CONFIDENCE < 0.8 ######

##### SHIFT IDS ######

# Inspect to see if any of the ranks are in the wrong spot due to missing taxonomic annotations - none are missing
rows_to_change_p <- !seqs_raw$phylum_rank %in% "phylum"
length(rows_to_change_p[rows_to_change_p==TRUE]) 

rows_to_change_c <- !seqs_raw$class_rank %in% "class"
length(rows_to_change_c[rows_to_change_c==TRUE]) 

rows_to_change_o <- !seqs_raw$order_rank %in% "order"
length(rows_to_change_o[rows_to_change_o==TRUE]) 

rows_to_change_f <- !seqs_raw$family_rank %in% "family"
length(rows_to_change_f[rows_to_change_f==TRUE])

rows_to_change_g <- !seqs_raw$genus_rank %in% "genus"
length(rows_to_change_g[rows_to_change_g==TRUE]) 

##### FILTER BY CONFIDENCE ######

## Remove any identifications that have a confidence < 0.8. Go from lowest taxonomic rank (species) to highest (phylum)

CO1v5_80 = seqs_raw %>% 
  mutate_at(c("species_name","species_confidence", "species_rank"), ~ifelse(species_confidence < 0.8,NA,.)) %>%
  mutate_at(c("genus_name","genus_confidence", "genus_rank"), ~ifelse(genus_confidence < 0.8,NA,.)) %>%
  mutate_at(c("family_name","family_confidence","family_rank"), ~ifelse(family_confidence < 0.8,NA,.)) %>%
  mutate_at(c("order_name","order_confidence","order_rank"), ~ifelse(order_confidence < 0.8,NA,.)) %>%
  mutate_at(c("class_name","class_confidence","class_rank"), ~ifelse(class_confidence < 0.8,NA,.)) %>%
  mutate_at(c("phylum_name","phylum_confidence","phylum_rank"), ~ifelse(phylum_confidence < 0.8,NA,.))%>%
  mutate_at(c("kingdom_name","kingdom_confidence","kingdom_rank"), ~ifelse(kingdom_confidence < 0.8, NA,.))


## And then, if there are blank names (e.g., not NA), change to NA
CO1v5_80 = CO1v5_80%>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

CO1v5_80_only_d = seqs_raw[seqs_raw$domain_confidence >= 0.8,]
CO1v5_80_only_k = seqs_raw[seqs_raw$kingdom_confidence >= 0.8,]
CO1v5_80_only_p = seqs_raw[seqs_raw$phylum_confidence >= 0.8,]
CO1v5_80_only_c = seqs_raw[seqs_raw$class_confidence >= 0.8,]
CO1v5_80_only_f = seqs_raw[seqs_raw$family_confidence >= 0.8,]
CO1v5_80_only_g = seqs_raw[seqs_raw$genus_confidence >= 0.8,]
CO1v5_80_only_s = seqs_raw[seqs_raw$species_confidence >= 0.8,]

##### MAKE GLOBAL ID AND CONFIDENCE COLUMN ######
# This will be the lowest taxonomic rank and confidence of an ASV ID in its own column 
# Makes it easier to filter later
# Also, distinguish lowest taxonomic rank 
CO1v5_80$confidence <- ifelse(!is.na(CO1v5_80$species_confidence), CO1v5_80$species_confidence,
                              ifelse(!is.na(CO1v5_80$genus_confidence), CO1v5_80$genus_confidence,
                                     ifelse(!is.na(CO1v5_80$family_confidence), CO1v5_80$family_confidence,
                                            ifelse(!is.na(CO1v5_80$order_confidence), CO1v5_80$order_confidence,
                                                   ifelse(!is.na(CO1v5_80$class_confidence), CO1v5_80$class_confidence,
                                                          ifelse(!is.na(CO1v5_80$phylum_confidence), CO1v5_80$phylum_confidence,
                                                                 ifelse(!is.na(CO1v5_80$kingdom_confidence), CO1v5_80$kingdom_confidence,
                                                                        ifelse(!is.na(CO1v5_80$domain_confidence), CO1v5_80$domain_confidence,  NA))))))))

CO1v5_80$id <- ifelse(!is.na(CO1v5_80$species_name), CO1v5_80$species_name,
                      ifelse(!is.na(CO1v5_80$genus_name), CO1v5_80$genus_name,
                             ifelse(!is.na(CO1v5_80$family_name), CO1v5_80$family_name,
                                    ifelse(!is.na(CO1v5_80$order_name), CO1v5_80$order_name,
                                           ifelse(!is.na(CO1v5_80$class_name), CO1v5_80$class_name,
                                                  ifelse(!is.na(CO1v5_80$phylum_name), CO1v5_80$phylum_name,
                                                         ifelse(!is.na(CO1v5_80$kingdom_name), CO1v5_80$kingdom_name,
                                                                ifelse(!is.na(CO1v5_80$domain_name), CO1v5_80$domain_name, NA))))))))

CO1v5_80$rank <- ifelse(!is.na(CO1v5_80$species_rank), CO1v5_80$species_rank,
                        ifelse(!is.na(CO1v5_80$genus_rank), CO1v5_80$genus_rank,
                               ifelse(!is.na(CO1v5_80$family_rank), CO1v5_80$family_rank,
                                      ifelse(!is.na(CO1v5_80$order_rank), CO1v5_80$order_rank,
                                             ifelse(!is.na(CO1v5_80$class_rank), CO1v5_80$class_rank,
                                                    ifelse(!is.na(CO1v5_80$phylum_rank), CO1v5_80$phylum_rank,
                                                           ifelse(!is.na(CO1v5_80$kingdom_rank), CO1v5_80$kingdom_rank,
                                                                  ifelse(!is.na(CO1v5_80$domain_rank), CO1v5_80$domain_rank,  NA))))))))

# Reorder these columns to the front
CO1v5_80 = CO1v5_80%>%
  dplyr::select(c("id","rank", "confidence"), everything())

# Set ASV as row name
rownames(CO1v5_80) <- CO1v5_80[,"name"]
CO1v5_80 = CO1v5_80[,-4]

writexl::write_xlsx(CO1v5_80, path = "CA_CO1_v5_80.xlsx")
