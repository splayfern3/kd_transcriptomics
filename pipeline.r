#load relevant libraries and illumina_loader file
library(tidyverse)
library(limma)
library(sva)
library(umap)
library(cluster)
library(msigdbr)
library(fgsea)
library(illuminaHumanv3.db)
library(illuminaHumanv4.db)
library(AnnotationDbi)
source("r_code/illumina_loader.R")

#load gse73461 to gse73463 along with the expression and pvalue files
gse73461 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt"
)
expr73461_lin <- gse73461$expr_lin
pval73461 <- gse73461$pval

expr73461_log2 <- log2(expr73461_lin + 1) #The paper used RSN normalization which does not do log2 transformation so we must do it ourself

#gse73462 and 73463 do not undergo this log2 transformation because they were already transformed when inputted into RSN
gse73462 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73462_GEOupload_Validation_HT12V3_Dataset_Normalised_Sept_15_n_147.txt"
)
expr73462_log2 <- gse73462$expr_lin 
pval73462 <- gse73462$pval


gse73463 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73463_GEOupload_Validation_HT12V4_Dataset_Normalised_Sept_15_n_233.txt"
)
expr73463_log2 <- gse73463$expr_lin
pval73463 <- gse73463$pval

#saving the rows of log2 normalized GSE data as probe IDs which will be mapped by AnnotationDbi 
probe_ids61 <- rownames(expr73461_log2) #uses illuminaHumanv4 chip
probe_ids62 <- rownames(expr73462_log2) #uses illuminaHumanv3 chip
probe_ids63 <- rownames(expr73463_log2) #uses illuminaHumanv4 chip

gplv4 <- read.delim("transcriptome_data/platforms/v4/GPL10558.txt", comment.char = "#", stringsAsFactors = FALSE)
gplv3 <- read.delim("transcriptome_data/platforms/v3/GPL6947.txt", comment.char = "#", stringsAsFactors = FALSE)
#map probe ids to gene name and symbol 
annov4 <- gplv4[c("Array_Address_Id", "Symbol", "Definition")]
annov3 <- gplv3[c("Array_Address_Id", "Symbol", "Definition")]



