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
library("WGCNA")
source("r_code/illumina_loader.R")
source("r_code/gene_collapse.R")
source("r_code/meta_filter.R")
source("r_code/process_study.R")
# TODO: add getGEO w/ conditionals so pipeline can be ran w/o assuming files already downloaded
# gse73463_meta <- getGEO("GSE73463", destdir = "transcriptome_data/uncompressed", getGPL = FALSE)




#load gse73461 to gse73463 along with the expression and pvalue files
gse73461 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt"
)
expr73461_lin <- gse73461$expr_lin
expr73461_lin[expr73461_lin < 0] <- 0
pval73461 <- gse73461$pval
#Remove columns filled with NA
na_cols61 <- which(colSums(is.na(expr73461_lin)) == nrow(expr73461_lin))
if(length(na_cols61) > 0) expr73461_lin <- expr73461_lin[, -na_cols61]
  
expr73461_log2 <- log2(expr73461_lin + 1) #The paper used RSN normalization which does not do log2 transformation so we must do it 
gse73461_meta <- getGEO("GSE73461", destdir = "transcriptome_data/uncompressed", getGPL = FALSE) # get meta

#gse73462 and 73463 do not undergo this log2 transformation because they were already transformed when inputted into RSN
gse73462 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73462_GEOupload_Validation_HT12V3_Dataset_Normalised_Sept_15_n_147.txt",
  floor_negatives_to_zero = FALSE
)
expr73462_log2 <- gse73462$expr_lin 
pval73462 <- gse73462$pval
gse73462_meta <- getGEO("GSE73462", destdir = "transcriptome_data/uncompressed", getGPL = FALSE) # get meta

gse73463 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73463_GEOupload_Validation_HT12V4_Dataset_Normalised_Sept_15_n_233.txt",
  floor_negatives_to_zero = FALSE
)
expr73463_log2 <- gse73463$expr_lin
pval73463 <- gse73463$pval
expr73463_log2[expr73463_log2 < 0] <- 0 #This file had artifacts of -24 which are biologically reelevent this removes them and sets them to 0
gse73463_meta <- getGEO("GSE73463", destdir = "transcriptome_data/uncompressed", getGPL = FALSE) # get meta

#CLEANING

# For GSE73461
colnames(expr73461_log2) <- gsub("_Detection_Pval|_Detection Pval", "", colnames(expr73461_log2))
colnames(pval73461)      <- gsub("_Detection_Pval|_Detection Pval", "", colnames(pval73461))
meta61_full <- pData(gse73461_meta[[1]])
keep61 <- apply(meta61_full, 1, function(row) any(grepl("Kawasaki|Control|KD", row, ignore.case = TRUE)))
res61_meta_clean <- meta61_full[keep61, ]
rownames(res61_meta_clean) <- res61_meta_clean$title

# For GSE73462
colnames(expr73462_log2) <- gsub("_Detection_Pval|_Detection Pval", "", colnames(expr73462_log2))
colnames(pval73462)      <- gsub("_Detection_Pval|_Detection Pval", "", colnames(pval73462))
meta62_full <- pData(gse73462_meta[[1]])
# Ensure rownames match title
rownames(meta62_full) <- meta62_full$title


# For GSE73463 (You actually already have this for 73463, which is great!)
colnames(expr73463_log2) <- gsub("_Detection_Pval|_Detection Pval", "", colnames(expr73463_log2))
colnames(pval73463)      <- gsub("_Detection_Pval|_Detection Pval", "", colnames(pval73463))
meta63_full <- pData(gse73463_meta[[1]])
# Ensure rownames match title
rownames(meta63_full) <- meta63_full$title
#filter out rows with high p values (i.e their intensity is comparable to the background)

# --- Process GSE73461 ---
res61 <- filter_by_metadata("GSE73461", expr73461_log2, pval73461, threshold = 0.5)
common61 <- intersect(colnames(expr73461_log2), rownames(res61_meta_clean))
res61$expr <- expr73461_log2[, common61]
res61$metadata <- res61_meta_clean[common61, ]
# Note: if filter_by_metadata returns pvals, align them too:
# res61$pval <- pval73461[, common61]

# --- Process GSE73462 ---
res62 <- filter_by_metadata("GSE73462", expr73462_log2, pval73462, threshold = 0.5)
common62 <- intersect(colnames(expr73462_log2), rownames(meta62_full))
res62$expr <- expr73462_log2[, common62]
res62$metadata <- meta62_full[common62, ]

# --- Process GSE73463 ---
res63 <- filter_by_metadata("GSE73463", expr73463_log2, pval73463, threshold = 0.5)
common63 <- intersect(colnames(expr73463_log2), rownames(meta63_full))
res63$expr <- expr73463_log2[, common63]
res63$metadata <- meta63_full[common63, ]


#saving the rows of log2 normalized GSE data as probe IDs which will be mapped by AnnotationDbi 
probe_ids61 <- rownames(res61) #uses illuminaHumanv4 chip
probe_ids62 <- rownames(res62) #uses illuminaHumanv3 chip
probe_ids63 <- rownames(res63) #uses illuminaHumanv4 chip

gplv4 <- read.delim("transcriptome_data/platforms/v4/GPL10558.txt", comment.char = "#", stringsAsFactors = FALSE)
gplv3 <- read.delim("transcriptome_data/platforms/v3/GPL6947.txt", comment.char = "#", stringsAsFactors = FALSE)
#map probe ids to gene name and symbol 
annov4 <- gplv4[, c("Array_Address_Id", "Entrez_Gene_ID", "Symbol")]
annov3 <- gplv3[, c("Array_Address_Id", "Entrez_Gene_ID", "Symbol")]
#rename cols
colnames(annov4) <- c("PROBEID", "ENTREZID", "SYMBOL")
colnames(annov3) <- c("PROBEID", "ENTREZID", "SYMBOL")

#filtered files are mapped to ENTREZ IDs 
mapped61_for_collapse <- annov4 %>%
  filter(PROBEID %in% rownames(res61$expr)) %>%
  filter(!is.na(ENTREZID) & ENTREZID != "")
mapped62_for_collapse <- annov3 %>%
  filter(PROBEID %in% rownames(res62$expr)) %>%
  filter(!is.na(ENTREZID) & ENTREZID != "")
mapped63_for_collapse <- annov4 %>%
  filter(PROBEID %in% rownames(res63$expr)) %>%
  filter(!is.na(ENTREZID) & ENTREZID != "")



# Process GSE73461
# We pass the filtered expression matrix and the mapping table
expr61_final <- process_study(
  expr_mat = res61$expr, 
  metadata = res61$metadata, 
  mapped_df = mapped61_for_collapse, 
  dataset_name = "GSE73461"
)

expr62_final <- process_study(
  expr_mat = res62$expr, 
  metadata = res62$metadata, 
  mapped_df = mapped62_for_collapse, 
  dataset_name = "GSE73462"
)

expr63_final <- process_study(
  expr_mat = res63$expr, 
  metadata = res63$metadata, 
  mapped_df = mapped63_for_collapse, 
  dataset_name = "GSE73463"
)