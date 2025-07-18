###______________________________________________________________________ 
###
### SMITHSONIAN TROPICAL RESEARCH INSTITUTE
### STRI-SENACYT INTERNSHIP, 2021 - 2022
### ROHR REEF RESILIENCE PROJECT
### Coral Reef eDNA MONITORING
###
### BEPE INTERNSHIP, 2024 - 2025
### LAGEM - LABORATÓRIO DE GENÔMICA E ECOLOGIA MARINHA
###
### by Matthieu Leray, edited by Viviane G. Ali S. and Isabela Buxbaum
###______________________________________________________________________



# Sections 
# (1) Set up
# (2) Sequencing Data Processing 
# (3) Inferring Amplification Sequence Variants with DADA 
# (4) Taxonomic Assignment
# (5) Create phyloseq object 
# (6) Remove Contaminants 
# (7) Exploratory Graphs
# (8) Clustering ASVs into 97% OTUs
# (9) OTU curation
# (10) New Phyloseq Object
# (11) More Exploratory Graphs



#_________________________________________________________________________
#
#### Set up ####
#_________________________________________________________________________


#activating and checking packages 
library(dada2)
library(tidyverse)
library(phyloseq)
library(decontam)
library(DECIPHER)
library(speedyseq)
library(metagMisc)
library(lulu)
library(pairwiseAdonis)

#make sure it's all downloaded
packageVersion('package name')



#_________________________________________________________________________
#
#### Sequencing Data Processing ####
#_________________________________________________________________________


#creating file path to data
path <- "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\ROHR42_fastq_raw"
list.files (path)

##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files (path, pattern = "_R1_001.trimmed.fastq")) #I edited it here: in the original script it was "_R1_001.fastq", now it is "_R1_001.trimmed.fastq", which means i'm working with the sequences that were processed by Cutadapt.
fnRs <- sort(list.files (path, pattern = "_R2_001.trimmed.fastq"))
sample.names <- sapply(strsplit(fnFs, "_"), '[', 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#plotting quality profiles
qprofile_rev <- print (plotQualityProfile(fnFs, aggregate = TRUE)       
                       + ggtitle ("Forward"))                           
qprofile_rev <- print (plotQualityProfile(fnRs, aggregate = TRUE)       
                       
                       + ggtitle ("Reverse"))



#placing filtered files in a new filtered subdirectory
filtFs <- file.path (path, "filtered", paste0 (sample.names, "_F_filt.fastq"))
filtRs <- file.path (path, "filtered", paste0 (sample.names, "_R_filt.fastq"))
names (filtFs) <- sample.names
names (filtRs) <- sample.names

#filtering and trimming, here truncation at 190 (Fwd) and 130 (Rev) bp,
#2expected errors max (N discarded automatically)



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(190, 130),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)  #use (multithread = FALSE) in Windows

print(out)




#My results:
#                                        reads.in reads.out    #it takes under consideration the PAIR of reads, and not only the R1, it analyzes R2 as well even though it doesn't appear in the table below.
#AI-23-B10_S76_L001_R1_001.trimmed.fastq    28541     24610
#AI-23-B11_S84_L001_R1_001.trimmed.fastq    61369     54957
#AI-23-B12_S92_L001_R1_001.trimmed.fastq   110338    100538
#AI-23-B16_S29_L001_R1_001.trimmed.fastq    22587     20521
#AI-23-B17_S38_L001_R1_001.trimmed.fastq    33756     30704
#AI-23-B18_S45_L001_R1_001.trimmed.fastq    33064     29742
#AI-23-B22_S77_L001_R1_001.trimmed.fastq    18444     13784
#AI-23-B23_S85_L001_R1_001.trimmed.fastq   151953    137693
#AI-23-B24_S93_L001_R1_001.trimmed.fastq   100183     90582
#AI-23-B28_S30_L001_R1_001.trimmed.fastq    15394     14120
#AI-23-B29_S37_L001_R1_001.trimmed.fastq    31248     29102
#AI-23-B30_S46_L001_R1_001.trimmed.fastq    29512     26962
#AI-23-B34_S78_L001_R1_001.trimmed.fastq    86564     79285
#AI-23-B35_S86_L001_R1_001.trimmed.fastq   108802    101247
#AI-23-B36_S7_L001_R1_001.trimmed.fastq     41366     37012
#AI-23-B4_S28_L001_R1_001.trimmed.fastq     53079     48257
#AI-23-B5_S36_L001_R1_001.trimmed.fastq     50786     45602
#AI-23-B6_S44_L001_R1_001.trimmed.fastq     33251     29758
#AI-23-T1_S4_L001_R1_001.trimmed.fastq      47074     42429
#AI-23-T13_S5_L001_R1_001.trimmed.fastq     13459     12396
#AI-23-T14_S13_L001_R1_001.trimmed.fastq    17988     16473
#AI-23-T15_S21_L001_R1_001.trimmed.fastq     9734      8546
#AI-23-T19_S53_L001_R1_001.trimmed.fastq    41630     37881
#AI-23-T2_S12_L001_R1_001.trimmed.fastq     42862     38203
#AI-23-T20_S61_L001_R1_001.trimmed.fastq    33852     29846
#AI-23-T21_S69_L001_R1_001.trimmed.fastq    59360     52801
#AI-23-T25_S6_L001_R1_001.trimmed.fastq      5796      2927
#AI-23-T26_S14_L001_R1_001.trimmed.fastq    23730     21367
#AI-23-T27_S22_L001_R1_001.trimmed.fastq    30907     28133
#AI-23-T3_S20_L001_R1_001.trimmed.fastq     47195     42238
#AI-23-T31_S54_L001_R1_001.trimmed.fastq    41124     37696
#AI-23-T32_S62_L001_R1_001.trimmed.fastq    34301     30172
#AI-23-T33_S70_L001_R1_001.trimmed.fastq    38374     31627
#AI-23-T7_S52_L001_R1_001.trimmed.fastq     32615     29123
#AI-23-T8_S60_L001_R1_001.trimmed.fastq     32274     28978
#AI-23-T9_S68_L001_R1_001.trimmed.fastq     33815     30397
#AV-23-B10_S73_L001_R1_001.trimmed.fastq    68943     48342
#AV-23-B11_S81_L001_R1_001.trimmed.fastq    80934     68397
#AV-23-B12_S89_L001_R1_001.trimmed.fastq    63187     55948
#AV-23-B16_S26_L001_R1_001.trimmed.fastq    27783     24821
#AV-23-B17_S34_L001_R1_001.trimmed.fastq    22600     20233
#AV-23-B18_S42_L001_R1_001.trimmed.fastq    54769     49755
#AV-23-B22_S74_L001_R1_001.trimmed.fastq    85376     74171
#AV-23-B23_S82_L001_R1_001.trimmed.fastq   105847     96442
#AV-23-B24_S90_L001_R1_001.trimmed.fastq    56583     52074
#AV-23-B28_S27_L001_R1_001.trimmed.fastq    36506     33379
#AV-23-B29_S35_L001_R1_001.trimmed.fastq    63277     56569
#AV-23-B30_S43_L001_R1_001.trimmed.fastq    41625     37152
#AV-23-B34_S75_L001_R1_001.trimmed.fastq    31356     27488
#AV-23-B35_S83_L001_R1_001.trimmed.fastq    49766     44100
#AV-23-B36_S91_L001_R1_001.trimmed.fastq    20069     17700
#AV-23-B4_S25_L001_R1_001.trimmed.fastq     30249     27387
#AV-23-B5_S33_L001_R1_001.trimmed.fastq     36276     32414
#AV-23-B6_S41_L001_R1_001.trimmed.fastq     30999     27023
#AV-23-T1_S1_L001_R1_001.trimmed.fastq      23877     21186
#AV-23-T13_S2_L001_R1_001.trimmed.fastq     17627     15268
#AV-23-T14_S10_L001_R1_001.trimmed.fastq    20909     18530
#AV-23-T15_S18_L001_R1_001.trimmed.fastq    26490     24516
#AV-23-T19_S50_L001_R1_001.trimmed.fastq    53275     47582
#AV-23-T2_S9_L001_R1_001.trimmed.fastq      15911     11581
#AV-23-T20_S58_L001_R1_001.trimmed.fastq       12         6
#AV-23-T21_S66_L001_R1_001.trimmed.fastq   108599     98430
#AV-23-T25_S3_L001_R1_001.trimmed.fastq     40479     37405
#AV-23-T26_S11_L001_R1_001.trimmed.fastq    24847     17743
#AV-23-T27_S19_L001_R1_001.trimmed.fastq    55205     49153
#AV-23-T3_S17_L001_R1_001.trimmed.fastq     29249     26334
#AV-23-T31_S51_L001_R1_001.trimmed.fastq    59417     53672
#AV-23-T32_S59_L001_R1_001.trimmed.fastq    30358     25573
#AV-23-T33_S67_L001_R1_001.trimmed.fastq    62223     55671
#AV-23-T7_S49_L001_R1_001.trimmed.fastq     26728     23942
#AV-23-T8_S57_L001_R1_001.trimmed.fastq        17        11
#AV-23-T9_S65_L001_R1_001.trimmed.fastq     42311     36315
#CPCR1_S94_L001_R1_001.trimmed.fastq           45        24
#CPCR2_S95_L001_R1_001.trimmed.fastq           31        25
#CPCR3_S96_L001_R1_001.trimmed.fastq           22        14
#QV-23-B10_S79_L001_R1_001.trimmed.fastq    18153     14677
#QV-23-B11_S87_L001_R1_001.trimmed.fastq    38076     34277
#QV-23-B12_S8_L001_R1_001.trimmed.fastq     38783     32163
#QV-23-B16_S40_L001_R1_001.trimmed.fastq    65028     59297
#QV-23-B17_S48_L001_R1_001.trimmed.fastq    35544     32440
#QV-23-B18_S56_L001_R1_001.trimmed.fastq    57629     53220
#QV-23-B22_S88_L001_R1_001.trimmed.fastq    42700     38933
#QV-23-B4_S31_L001_R1_001.trimmed.fastq        21        12
#QV-23-B5_S39_L001_R1_001.trimmed.fastq     62936     57933
#QV-23-B6_S47_L001_R1_001.trimmed.fastq     28309     25922
#QV-23-T1_S15_L001_R1_001.trimmed.fastq        21        13
#QV-23-T13_S16_L001_R1_001.trimmed.fastq    42375     38294
#QV-23-T14_S24_L001_R1_001.trimmed.fastq    51926     47922
#QV-23-T15_S32_L001_R1_001.trimmed.fastq    65921     60044
#QV-23-T19_S64_L001_R1_001.trimmed.fastq    48356     44092
#QV-23-T20_S72_L001_R1_001.trimmed.fastq    59461     53613
#QV-23-T21_S80_L001_R1_001.trimmed.fastq    55410     50090
#QV-23-T3_S23_L001_R1_001.trimmed.fastq        25        18
#QV-23-T7_S55_L001_R1_001.trimmed.fastq     39427     36110
#QV-23-T8_S63_L001_R1_001.trimmed.fastq     44703     39586
#QV-23-T9_S71_L001_R1_001.trimmed.fastq    117300    106459



#to just consider files that pass the filtering 
exists <- file.exists (filtFs) & file.exists (filtRs)
filtFs <- filtFs [exists]
filtRs <- filtRs [exists]

# Learning error rates 
errF <- learnErrors (filtFs, multithread = TRUE)
errR <- learnErrors (filtRs, multithread = TRUE)


#My results:
# > errF <- learnErrors (filtFs, multithread = TRUE)
# 103807070 total bases in 546353 reads from 11 samples will be used for learning the error rates.
# > errR <- learnErrors (filtRs, multithread = TRUE)
# 102811670 total bases in 790859 reads from 15 samples will be used for learning the error rates.


# Plotting errors
plotErrors (errF, nominal = TRUE)
plotErrors (errR, nominal = TRUE)

# Dereplicating reads
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
derepFs <- derepFastq (filtFs)
names (derepFs) <- sam.names
derepRs <- derepFastq (filtRs)
names (derepRs) <- sam.names

# Save copies of dereplicated objects 
saveRDS(derepFs, file ="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\derepFs.rds")
saveRDS(derepRs, file ="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\derepRs.rds")


#____________________________________________________________________________
#
#### Infering Amplification Sequence Variants with dada ####
#____________________________________________________________________________

dadaFs <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
dadaFs [[1]]
dadaFs [[5]]
dadaFs [[96]] #You can choose any number among the number of samples you have, preferably the last sample.


#My results:
# > dadaFs <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
# 96 samples were pooled: 3617200 reads in 484680 unique sequences.
# > dadaFs [[1]]
# dada-class: object describing DADA2 denoising results
# 306 sequence variants were inferred from 4817 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
# > dadaFs [[5]]
# dada-class: object describing DADA2 denoising results
# 456 sequence variants were inferred from 6595 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
# > dadaFs [[96]] 
# dada-class: object describing DADA2 denoising results
# 446 sequence variants were inferred from 17501 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


dadaRs <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
dadaRs [[1]]
dadaRs [[2]]
dadaRs [[96]] #You can choose any number among the number of samples you have, preferably the last sample.


#My results:
# > dadaRs <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
# 96 samples were pooled: 3617200 reads in 394316 unique sequences.
# > dadaRs [[1]]
# dada-class: object describing DADA2 denoising results
# 285 sequence variants were inferred from 3839 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
# > dadaRs [[2]]
# dada-class: object describing DADA2 denoising results
# 531 sequence variants were inferred from 10994 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
# > dadaRs [[75]] #Here it said 100, but i changed it to 75 bcs its the number of samples i have
# dada-class: object describing DADA2 denoising results
# 414 sequence variants were inferred from 14158 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


# Save copies of dada objects
saveRDS(dadaFs, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\dadaFs_ITS2.rds")
saveRDS(dadaRs, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\dadaRs_ITS2.rds")


##Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtab1 <- makeSequenceTable(mergers)
dim(seqtab1)

#My results:
# >dim(seqtab1)
# [1]   96 2899
# This means in 96 samples we have 2899 different ASVs

# Inspect the merged sequences from the data frame of the first sample (and the
# 6th sample).
head(mergers[[1]])
head(mergers[[6]])


#My results
# > head(mergers[[1]])
# sequence
# 3  AGTTGGGTGCTTTCAGCACCGGACTATGCCTGGTTCAGACCCGCTTACAAACCAATCACACACATCTAATGTGTTGTGGACAATGGCAGTGACTCAACAGGCTGTTGGTCTCTGCTCAAGCGCAAGCCATCGTCTCGGTTAACCGGCTTCGGCCGGGGAGCCTCTGAGGCGTGACTTAGGTGATACCATCTTCGTCGCACCACCTCTAGGGAACTTTCAAAGACTAAACATCTTGCTAAAACATGGGCTGAGCTCAGGTAAGAGGACCCGC
# 7                                                GCTTGGGATGAACAGTCCCAAGCATGTCTGTTTGAGGGTCATTTTTATCTCACCCACCCTTGTGGTGTGTGGAGTGTGGCTGTCCTTTTGTAGAGAGAGGACGGCTCAAGATCAGGCCGACGAGTCTCTGGCAAGCGAGGCGACGCGAAGAGGAATCTTACGGCGCGGAGCGAGCACGGACGACGCAACTTGAATTATGACCTCAAATCAGGCAAGACTACCCGC
# 8  AGTTGGGTGCTTTCAGCACCGGACTATGCCTGGTTCAGACCCGCTTACAAACCGATCACACACATCTAATGTGTCGTGGACAATGGCAGTGACTCAACAGGCTGTTGGTCTCTGCTCAAGTGCAAGCCATCGTCTCGGTTAACCGGCTTCGGCCGGGGAGCCTCTGAGGCGTGACTTAGGTGATACCATCTTCGTCGCACCACCTCTAGGGAACTTTCAAAGACTAAACATCTTGCTAAAACATGGGCTGAGCTCAGGTAAGAGGACCCGC
# 12                                             GCTTGGGATGAACAGTCCCAAGCATGTCTGTTTGAGGGTCATTTTTATCTCACCCACCCCTGTGGTGTGTGGAGTGTGGCTGTCCTCTTTCTTAGAGAGAGGACGGCTCAAGATCAGGCCGACGAGTCTCTGGCAAGCGAGGCGACGCGAAGAGGAATCTTACGGCGCGGAGCGAGCACGGACGACGCAACTTGAATTATGACCTCAAATCAGGCAAGACTACCCGC
# 14        ACTGAGTTTCGGCTCAGTATGCCTGGGTCAGTCTCACTTACAACACCTCAACCATTCATTTGGTTGGACAATGGGAGTCAGGGGCGAACCGCTGCATGTCTCCTAAAAGAAAGGTTAACGCATCATCGTTGTGTCACCTTGTGTGGCGCGGCGCCGTGCCGCTGTCGGCTGGTTAACACCTCGACAGGGGAGGGCCCGTGTCCGTGACCATGTTGATTCGTCTTCTCCTTTCCCTTCGAGCTGATCTCAGGTAAGAACACCCGC
# 15                                           GCTTGGGATGAACAGTCCCAAGCATGTCTGTTTGAGGGTCATTTTTATCTCACCCACCCTCGTGGTGTGTGGAGTGTGGCTGTCCTCTTTCTTCGAGAGAGAGGACGGCTCAAGATCAGGCCGACGAGTCTCTGGCAAGCGAGGCGACGCGAAGAGGAATCTTACGGCGCGGAGCGAGCACGGACGACGCAACTTGAATTATGACCTCAAATCAGGCAAGACTACCCGC
# abundance forward reverse nmatch nmismatch nindel prefer accept
# 3       2551       2       1     49         0      0      2   TRUE
# 7        705       7       2     95         0      0      2   TRUE
# 8        579     138       1     49         0      0      2   TRUE
# 12       244      28       2     93         0      0      2   TRUE
# 14       203      57      46     56         0      0      2   TRUE
# 15       196      30       2     91         0      0      2   TRUE

# > head(mergers[[6]])
# 2         GCTCCCGGGGTCTCCTCGGGAGCAGGTCTGTCCGAGCGTCTGCCAACCATGAGTTCTCGCCGAGCGCGAGGAGATTGGGGAGCTGCGGCATCGTCGGGTACGACCGCATGCCGCATTCCCTGAAGATTTGACGACCGGTGTGGTCGCCCGCCGTACGAGACGGTCGGTGACAACGGTGCACAGTGTTGTGCCGCCGTGCTCGTGGGGTGTTGGCGGAGCGCCGGTGGAGTGCGATCAGGTCGTGCTTTTGCCAATCTGGACCTCGGATCAGGCCAGGCTACCCGC
# 3                                                                     GCTTGGGATGAACAGTCCCAAGCATGTCTGTTTGAGGGTCATTTTTATCTCACCCACCCTTGTGGTGTGTGGAGTGTGGCTGTCCTTTTGTAGAGAGAGGACGGCTCAAGATCAGGCCGACGAGTCTCTGGCAAGCGAGGCGACGCGAAGAGGAATCTTACGGCGCGGAGCGAGCACGGACGACGCAACTTGAATTATGACCTCAAATCAGGCAAGACTACCCGC
# 4                                    GCTCTCTGGATATCCAGGGAGCATGCCTGTCTGAGCGTCGTTTAACAATCTACACACAATTTGTGTGGAAAGATCATTCATGCTTCAATGTGTGTCGATCTATGGATTTGAACCGTGTGGTACTCCTTTGTGTACTGCTTGGAAGGCTCAACTTGCTTCTGTCCTTCTCTGCAGTGTCATGATGTGTGTATAACGTGTTCTCTTTCGACCTAGTCGAAAATTCTTATTTTTGACCTCAGATCAGACAAGACTACCCGC
# 5 GCCCTCTACCTCGTGTCGAGGGCATGCCTGGATCAGTCTTACTTTTACTCTCCCAGCCACTTCGGTGATGGTTGAGAAGTTGGGCGTCACGCAATCGGACTATTCGCGTGTCGCTTTAAGGTTATAAGGTCTACTGACTGCTAGTGGTTCGGTCCCTCGAGGGCAACCTCGACGGGACACCCAAGAGCGCAGTTTCGGTTGGTCTCTGACTGACTGATACCGAGCTCAAGGTGAACTCTGGTTGGAAGTGTTCAACTTGTAACTCCAAGCTGATCTCAGGTAAGGACACCCGC
# 6                         GCTTCACGCTTGCGTGGAGCATGCCTGGGTCAGTCTAATTTACAACACCTCAACCATTCTATGGTTGGATAATGGAAGTCACGCACCGGACCGGTCACGTGTCTTCTCAATGAAAGGTCAAACGCATCACTTCAGTGTCACCTCGTGTGATTCAGAGTCATGCCGTTGTCGGTTGGGTAAAACCTCGACATACGGAAGACACTAATCCGTGACCGCTGTGATAGTAACTCTCCTTTCCTTTCAAGCTGATCTCAGGTAAGGACACCCGC
# 7                                    GCTCTCTGGATATCCAGGGAGCATGCCTGTCTGAGCGTCGTTTAACAATCTACACACAATTTGTGTGGAAAGATCATTCATGCTTCAATGTGTGTCGATCTATGGATTTGAACCGTGTGGTACTCCTTTGTGTACTGCTTGGAAGGCTCAACTTGCTTCTGTCCTTCTCTGCAGTGTCATGATGTGTGTATAACGTGTTCTCTTTCGACCTAGTCGAAAATTCTTACTTTTGACCTCAGATCAGACAAGACTACCCGC
# abundance forward reverse nmatch nmismatch nindel prefer accept
# 2      1143       5       5     35         0      0      2   TRUE
# 3       716       6       2     95         0      0      2   TRUE
# 4       706      12      15     62         0      0      1   TRUE
# 5       606      13      18     27         0      0      1   TRUE
# 6       461      42      27     51         0      0      2   TRUE
# 7       318      12      26     62         0      0      1   TRUE


table (nchar(getSequences(seqtab1)))


#My results:
# > table (nchar(getSequences(seqtab1)))
# 190 191 193 195 196 197 198 199 202 203 204 205 206 207 208 209 210 211 213 214 215 217 218 220 221 222 223 224 225 227 228 229 230 231 232 233 234 
#  54  11   2   1   3   2   2   1   1   1   1   2   4   5   8   2   3   1   3   5   1   1   1   1   8   4   6   2  19  12   5 220   2   3   6   9   4 
# 235 236 237 238 239 240 241 242 243 244 245 246 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 
#   3   1   1   3   1   7   4   3   8   7   7  19  10   7   4   8  61  28  17   7   3   7  33  16  13  40  35  32  39 111  26  51  58 102  28  99 160 
# 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 
#  59  55  94  16   7  57  19  36  55  33  79  66 116  29  22  22  21  56  31  34  60  93  77  26   8   6  10  14  13  21  14  37  40  20  34  44

# This shows the size of the sequences we have and beneath it is the number of sequences that have this size


# Saving files to use in the next part of the workflow
saveRDS(seqtab1, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab1.rds")

# Save copies of seqtab (sequence table)
saveRDS(seqtab1, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab1.rds")
write.csv(seqtab1, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab+Chim.csv")


# Identifying and removing chimeras
seqtab1.nochim <- removeBimeraDenovo(seqtab1, 
                                     method="consensus", 
                                     multithread=TRUE)
dim (seqtab1.nochim)


#My results:
# > dim (seqtab1.nochim)
# [1]   96 1756
# This means in 96 samples we have 2899 different ASVs


# Saving table without chimeras and downloading as .csv file

saveRDS (seqtab1.nochim, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab_nochim.rds")
write.csv (seqtab1.nochim, file = "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\nochim_seq.csv")

# Generate .txt document with the changes generated in each step

getN <- function(x) sum(getUniques(x))
track <-    cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                  sapply(mergers, getN), rowSums(seqtab1.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sam.names
write.table(track, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\tracked_its2.txt")
saveRDS(seqtab1.nochim, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab1.nochim.rds")


#____________________________________________________________________________
#
#### Taxonomic Assignment ####
#____________________________________________________________________________


library(devtools)
library(tidyverse)
library(insect)

## Saving no-chimera csv as a fasta file. 
# Before this step you should:
# Create from scratch a new CSV file with two columns that contains the ASV IDs and their respective sequences (I called it reference_plato1_completo")
# You can do it by copying and pasting the sequences from nochim_seq.csv using the function "TRANSPOSE" in Google Sheet


setwd("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA")

source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
TabularToFasta("reference_plato1_completo.csv")
# Sometimes, when creating the fasta file it excludes the first line (ASV1), so pay attention to that in case you need to manually add it to the new fasta file
# I'm pretty sure that happens when you don't name the columns. So make sure to use the first line to name "ASVS, Sequences"


##Assign taxonomy with insect package and provided databases as reference: 
# (https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html)
# Download the classifier.rds from the marker you want and save it in your computer

# get sequences from table column names
x <-char2dna(colnames(seqtab1.nochim))
# name the sequences sequentially
names(x) <- paste0("ASVA", seq_along(x))  #I'm adding A, B, C, D... at the end of ASV so the ASV IDs from this analysis don't get mixed with other analyses when i put it all together in a single phyloseq obj from all 4 plates
# optionally remove column names that can flood the console when printed
# exchange sequences for their respective ASV name
colnames(seqtab1.nochim) <- names(x)
# save new table with modified names
write.csv(seqtab1.nochim, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\nochim_seq_new_names.csv", row.names = TRUE)

packageVersion("insect")
readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\seqtab1.nochim.rds")
setwd("C:/Users/isabu/Documents/Isabela/Oceanografia/LerayLab/Data_analyses/ITS2_Isabela") #this is where you saved your classifier.rds

classifier <- readRDS("classifier.rds")
x <- readFASTA("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\dna_fasta.fasta")
out <- classify(x, classifier, threshold = 0.8)

saveRDS(out, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\taxonomic_assignment.rds")
write.csv(out, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\taxonomic_assignment.csv")
View(out)


#________________________________________________________________________
#
#### Create phyloseq Object ####
#________________________________________________________________________

#Creating a phyloseq object - charge neccesary libraries
library(phyloseq)
library(tibble)
library(dplyr)
library(phyloseq)
library(Biostrings)

#At the end of this script, I found trouble concluding the taxonomic assignment because the original script did not include the csv with the reference sequences in the phyloseq obj.
#So I added a few more steps in order to create another phyloseq containing the ref seqs, which is the one we will be using afterwards to continue with the analysis

#Additionally, the csv files automatically generated previously by dada2 (ASVtable - nochim_seq_new_names.csv) and insect (TAXtable - taxonomic_assignment) are generated with the wrong format.
#So make sure to create new csv files, containing the same information, naming the necessary columns with 'SampleID', 'ASVS' and excluding any extra (and unnecessary) column/line that you might come across.

#Also, create a new csv file containing all the metadata from your data (it's the SAMtable)

#Creating phyloseq object - reading the 3 csv tables: 
#ASVtable from dada2
#TAXtable from insect
#SAMtable - a new table with the variable information of the samples
asvtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\nochim_seq_new_names_useful.csv")
taxtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\taxonomic_assignment_useful.csv")
samtable <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\SAMtable_plato4.csv")

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
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\ROHR42.1_obj.rds")

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
seqfile <- "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\reference_plato4.csv"
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
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\full_ROHR42.1_obj.rds")



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
          file = "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_contamdf.prev.csv")

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
write.csv(df.pa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\df.pa.csv", row.names = FALSE)

#Save work to this point
save(ROHR42.1_obj, contamdf.prev, file = "my_important_objects.RData")
save.image("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\phyloseq-decontam.RData")


# Create a new phyloseq object without contaminants identified in contamdf.prev
# The original script removed SEQUENCES that were KNOWN as contaminants, since we don't know that explicitly up to this point, I did it using the information we've got so far.
# So I did it but automatically removing what the script interpreted as "TRUE" to contaminants in the previous steps, not by manually inserting the sequences


contamdf.prev <- isContaminant(ROHR42.1_obj, method="prevalence", neg="is.neg")

#Extract the names od the ASVs determined as contaminants

badTaxa <- rownames(contamdf.prev)[which(contamdf.prev$contaminant)]
length(badTaxa)         #shows how many ASVs are considered contaminantes
head(badTaxa, n = 5)    #shows contaminants

#Create new phyloseq without the contaminants
goodTaxa <- setdiff(taxa_names(ROHR42.1_obj), badTaxa)
ROHR42.1_obj <- prune_taxa(goodTaxa, ROHR42.1_obj)

#Now remove all control samples 
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR10")
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR11")
ROHR42.1_obj = subset_samples(ROHR42.1_obj, sample_names(ROHR42.1_obj) !="CPCR12")

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
##otu_table()   OTU Table:          [ 183 taxa and 5 samples ]:
##sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
##tax_table()   Taxonomy Table:     [ 183 taxa by 11 taxonomic ranks ]:
##refseq()      DNAStringSet:       [ 183 reference sequences ]
##taxa are columns


#Save phyloseq object
saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

#Reading the phyloseq object
ROHR42.1 <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

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

saveRDS(ROHR42.1_obj, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")

loaded_obj <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")
head(refseq(loaded_obj))

loaded_obj <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")
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
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Location~Site)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Location~Depth)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Location~Site~Depth)
plot_bar(na.rh1, x = "order", fill = "order", facet_grid = ~Location~Season~Depth)

#Prepare for ordination plots and do plots
na.rh1 <- prune_samples(sample_sums(na.rh1) > 0, na.rh1)
rh1.ord <- ordinate(na.rh1, "NMDS", "bray")
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Location", 
                title = "Samples by Locality") + geom_point(size = 2)
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Season", 
                title = "Samples by Season") + geom_point(size = 2)
plot_ordination(na.rh1, rh1.ord, type = "samples", color = "Depth", 
                title = "Samples by profundidade") + geom_point(size = 2)
rh1_ord2 <- ordinate(na.rh1, "PCoA", "bray")
plot_ordination(na.rh1, rh1_ord2, type = "samples", color = "Location", 
                title = "PCoA for Samples by Locality") + geom_point(size = 2)


#normalizing number of reads for relative abundance
total = median(sample_sums(na.rh1))
standf = function(x, t=total) round(t * (x/sum(x)))
na.rh1 = transform_sample_counts(na.rh1, standf)
na.rh1


#My results:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 7 taxa and 5 samples ]:
#   sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 7 taxa by 8 taxonomic ranks ]:
#   refseq()      DNAStringSet:       [ 7 reference sequences ]
# taxa are columns


#drawing plots with normalized data
plot_bar(na.rh1, fill = "order")
plot_bar(na.rh1, x = "Location", fill = "order")
plot_bar(na.rh1, x = "Sample", fill = "order") +
  facet_wrap(~Location)


#merging samples by locality for ploting
sample_data(na.rh1)$Site <- factor(sample_data(na.rh1)$Location)
sum(is.na(sample_data(na.rh1)$Location))
na.rh1 <- subset_samples(na.rh1, !is.na(Location))
na.rh1_locality <- merge_samples(na.rh1, "Location")
plot_bar(na.rh1_locality, fill = "order") 

plot_bar(na.rh1, x = "Site", fill = "order") +
  facet_grid(. ~ Location)  # ou facet_grid(Location ~ .)


#plotting heatmap for abundance by locality
otu_table(na.rh1) <- otu_table(otu_table(na.rh1) + 1, taxa_are_rows = TRUE)
(p <- plot_heatmap(na.rh1, "NMDS", "bray", "Location", "order"))

p <- plot_heatmap(na.rh1, method = "NMDS", distance = "bray", sample.label = "Location", taxa.label = "order")
print(p)  #Força a exibição


#Filtering data to consider taxa present at least at 20%
rh1_abund <- filter_taxa(na.rh1, function(x) sum(x > total*0.20) > 0, TRUE)
rh1_abund


#My results:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 4 taxa and 5 samples ]:
#   sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 4 taxa by 8 taxonomic ranks ]:
#   refseq()      DNAStringSet:       [ 4 reference sequences ]
# taxa are columns


#run heatmap for filtered data
(q <- plot_heatmap(rh1_abund, "NMDS", "bray", "Location", "order"))
q <- plot_heatmap(na.rh1, method = "NMDS", distance = "bray", sample.label = "Location", taxa.label = "order")
print(q)

#ordination plots with normalized data
#generating object for nmds with bray curtis and ploting
normalize.ord <- ordinate(na.rh1, "NMDS", "bray")
plot_ordination(na.rh1, normalize.ord, type = "samples", color = "Location", 
                title = "Samples by Locality") + geom_point(size = 3)

#generating object for PCoA with bray curtis and ploting
normal.pcoa <- ordinate(na.rh1, method = "PCoA", distance = "bray")
plot_ordination(na.rh1, normal.pcoa, type = "samples", color = "Location", 
                title = "PCoA for Samples by Locality") +
  geom_point(size = 4)


#phyloseq graphs...




#________________________________________________________________________
#
#### Clusterinng ASVs into 97% OTUs ####
#________________________________________________________________________



setwd("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu")

nproc <- 8
its2 <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\ROHR42.1_final_obj.rds")
sequences <- taxa_names(its2)
sample_names <- sample_names(its2)

#Verifies if refseq is present
if (is.null(refseq(its2))) {
  stop("O slot 'refseq' está vazio. Verifique se o objeto salvo contém as sequências.")
}


#Installs decipher 
dna <- refseq(its2) 
set.seed(123) #initialize the random number generator
clusters <- DECIPHER::Clusterize(dna, method = "overlap", 
                                 cutoff = 0.03, #97% 
                                 penalizeGapLetterMatches = NA, 
                                 includeTerminalGaps = TRUE, 
                                 processors = nproc)
set.seed(NULL)


#penalize gap-to-letter mismatches once per insertion or deletion, 
#which treats runs of gaps (i.e., indels) as equivalent to a single mismatch
#the calculation of distance will use the entire (global) alignment
#Note: function "merge_taxa_vec" was in the package mikemc/speedyseq 
#NOTE VIVIANE: 99% OTUS gave 6004 ASV's | 97% OTUS gave 2408 ASV's --> I'm choosing 97% otus for the analysis. 
#Given that diversity obtained through both approaches is basically the same. 


#Merge clusters with taxonomy
its2.otu <- merge_taxa_vec(
  its2,
  group = clusters$cluster,
  tax_adjust = 0)


#Checking objects 

its2


#My results:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 183 taxa and 5 samples ]:
# sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
# tax_table()   Taxonomy Table:     [ 183 taxa by 11 taxonomic ranks ]:
# refseq()      DNAStringSet:       [ 183 reference sequences ]
# taxa are columns


its2.otu #97%

#My results:
#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 110 taxa and 5 samples ]:
#sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
#tax_table()   Taxonomy Table:     [ 110 taxa by 11 taxonomic ranks ]:
#refseq()      DNAStringSet:       [ 110 reference sequences ]
#taxa are columns


# Take a look at files
tax_dada_tax  = tax_table(its2.otu)
write.csv(tax_dada_tax, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\otu_tax.csv")

tax_dada_otu  = otu_table(its2.otu, taxa_are_rows = FALSE)
write.csv(tax_dada_otu, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\otu_matrix.csv")

#Extra step to generate a csv containing the reference sequences from the phyloseq obj
#Verifies if there's a refseq
refseqs_clustered <- refseq(its2.otu)

#Generates a dataframe with the ASV names (in this case it's the OTU names) and their respective sequences
asv_seq_df <- data.frame(
  ASVA = names(refseqs_clustered),  #I added a C at the end of ASV so these ASVs don't get mixed with the ASVs from other phyloseq obj when I put them all together
  Sequence = as.character(refseqs_clustered)
)

#Create a CSV
write.csv(asv_seq_df, file = "otu_sequences_97percent.csv", row.names = FALSE)



#________________________________________________________________________
#
#### OTU curation ####
#________________________________________________________________________


# File preparation
# A - OTU table with samples as columns and OTUs as rows
tax_dada_otu  = t(tax_dada_otu)
tax_dada_otu.df  = phyloseq_to_df(tax_dada_otu, addtax = F, addtot = F, addmaxrank = F,
                                  sorting = "abundance")  #with package metagMisc
tax_dada_otu.df  <- data.frame(tax_dada_otu.df , row.names = 1)

# B - Fasta file to prepare the match list
# from .csv file created in the previous topic
fa.table = read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\otu_sequences_97percent.csv")
fa = character(2 * nrow(fa.table))
fa[c(TRUE, FALSE)] = paste0(">", fa.table$ASV)
fa[c(FALSE, TRUE)] = as.character(fa.table$Sequence)

writeLines(fa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\otus97.fasta")

# C - Create match list with vsearch
# Place input file "otus97.fasta" in "/Applications/vsearch/bin" folder of vsearch and run:
# vsearch --usearch_global otus97.fasta --db otus97.fasta --self --id .84 --iddef 1 --userout match_list97.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

matchlist_name = read.table("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\vsearch-2.30.0-win-x86_64\\bin\\match_list.txt")
names(matchlist_name)[names(matchlist_name) == "V1"] <- "OTUid"
names(matchlist_name)[names(matchlist_name) == "V2"] <- "hit"
names(matchlist_name)[names(matchlist_name) == "V3"] <- "match"
matchlist_name$OTUid <- as.character(matchlist_name$OTUid)
matchlist_name$hit <- as.character(matchlist_name$hit)

# Run OTU curation
curated_result <- lulu(tax_dada_otu.df, matchlist_name)

# Curated OTU table
curated_table = curated_result$curated_table
write.csv(curated_table, file="C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\curated_table97.csv")

# Curation Results 
# 99 97% OTUS were obtained. 

##Prepare fasta file for taxonomic assignment using the curated table. 
#In order to do that, create a sequence list by saving first column of "curated_table97.csv" as a new file. Name the column "ASV(letter)" and beneath that put the ASV IDs.
#I called my csv file "curated97_seq.csv". 
#Since my curated_table97.csv contains only the ASV IDs and not their sequences, I added a new part to this code that makes the association of the ASV ID with their correspondent sequence so we can then create the fasta file.

#Read both csv files
asv_names <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\curated97_seq.csv", stringsAsFactors = FALSE)
asv_seqs <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\otu_sequences_97percent.csv", stringsAsFactors = FALSE)

#Check if both have "ASVA" as column names
names(asv_names)
names(asv_seqs)

#Merge them using the column ASVA
merged_df <- merge(asv_names, asv_seqs, by = "ASVA")

#Save new file
write.csv(merged_df, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\asv_com_sequencias.csv", row.names = FALSE)

##

curated_seqs = read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\asv_com_sequencias.csv")
fa = character(2 * nrow(curated_seqs))
fa[c(TRUE, FALSE)] = paste0(">", curated_seqs$ASV)
fa[c(FALSE, TRUE)] = as.character(curated_seqs$Sequence)
writeLines(fa, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\curated_97OTUseqs_its2.fasta")




#________________________________________________________________________
#
#### Repeat taxonomic assignment with insect ####
#________________________________________________________________________

# --> Done by Helio Quintero at NAOS Computers. 
# --> The classifier used was provided by insect github
# --> https://github.com/shaunpwilkinson/insect


x <- readFASTA("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\curated_97OTUseqs_its2.fasta")
classifier <- readRDS("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\classifier.rds") 
out <- classify(x, classifier, ping = 0.98, cores = 8)

write.csv(out, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\otus97_insect.csv")



#_________________________________________________________________________
#
#### New Phyloseq ####
#_________________________________________________________________________

#In order to do this, again, we need to correct the format of the files otus97_insect.csv and curated_table97.csv created in previous steps.
#So make sure to create new csv files, containing the same information, naming the necessary columns with 'SampleID', 'ASVS' and excluding any extra (and unnecessary) column/line that you might come across.
#As well as correcting the format of the otumat97.csv that will be generated and create a new csv containing the metadata.

# Transpose otu table
otu_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\curated_table97_useful.csv")
asv_mat <- t(otu_mat)
write.csv(asv_mat, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\otumat97.csv")

# Read tables
otu_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\otumat97_useful.csv")
tax_mat <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\otus97_insect_useful.csv")
sample_df <- read.csv("C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\SAMtable_otu_plato4.csv")

# Define rownames
colnames(otu_mat)

otu_mat <- otu_mat %>% 
  tibble::column_to_rownames("SampleID") 
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("ASVS") 
sample_df <- sample_df %>%
  tibble::column_to_rownames("SampleID") 

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
plato1_completo_otu <- phyloseq(asv, tax, sample)
plato1_completo_otu

#My results:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 103 taxa and 5 samples ]:
#   sample_data() Sample Data:        [ 5 samples by 8 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 103 taxa by 11 taxonomic ranks ]:
#   taxa are rows

head(taxa_names(asv))
head(taxa_names(tax))

saveRDS(plato1_completo_otu, "C:\\Users\\isabu\\Documents\\Isabela\\Oceanografia\\LerayLab\\Data_analyses\\ITS2_Isabela\\Fastq_ROHR42\\plato_1_completo\\Rscript\\1_SDP\\2_IASVD\\3_TA\\4_CPO\\5_RC\\6_EG\\7_clustering_asv_into_otu\\8_otu_curation\\9_NP\\plato4_otu.rds")


#Verify general structure

summary(plato1_completo_otu)


#Verify richness/diversity
#Separate the phyloseq obj by location
buzios_only <- subset_samples(plato1_completo_otu, Location == "Ilha de Búzios")
laje_only <- subset_samples(plato1_completo_otu, Location == "Laje de Santos")

#Plot richness/diversity to each island
plot_richness(laje_only, x = "Site", measures = c("Shannon", "Observed")) +
  ggtitle("Laje de Santos")

plot_richness(buzios_only, x = "Site", measures = c("Shannon", "Observed")) +
  ggtitle("Queimada Grande")


#Checar composição por filo, etc...
# Filtrar Alcatrazes
buzios_filtrado <- buzios_only %>%
  subset_taxa(!phylum %in% c("Ctenophora", "Nematoda") & !order %in% c("", "Philasterida"))

# Filtrar Laje de Santos
laje_filtrado <- laje_only %>%
  subset_taxa(!phylum %in% c("Ctenophora", "Nematoda") & !order %in% c("", "Philasterida"))

# Plotar barras por filo
plot_bar(buzios_filtrado, fill = "order") +
  ggtitle("Composição por ordem - Ilha de Búzios")

plot_bar(laje_filtrado, fill = "order") +
  ggtitle("Composição por ordem - Laje de Santos")



#Verificar riqueza/diversidade

plot_richness(buzios_filtrado, x = "Site", measures = c("Shannon", "Observed")) #Substitua "SEU_GRUPO" pelo nome de uma variável categórica em sample_data, como Site, Depth, Treatment, etc.


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
#saveRDS(rohr_out, 
"C:/Users/Viviane/Documents/2024/Benthic_eDNA/rohr_out.rds")


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


