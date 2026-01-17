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
library(ggplot2)
library(harmony)
source("r_code/illumina_loader.R")
source("r_code/gene_collapse.R")
source("r_code/meta_filter.R")
source("r_code/process_study.R")
source("r_code/clean_dx.R")

install.packages(c("devtools", "clv", "fields", "matrixStats", "data.table", "cluster", "clue", "circlize", "gdata"))

# Install the Bioconductor dependency
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

# Finally, install the RCC package itself
devtools::install_github("MSCTR/RecursiveConsensusClustering"

library(RecursiveConsensusClustering)
                         
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

common_genes <- intersect(rownames(expr61_final), 
                          intersect(rownames(expr62_final), rownames(expr63_final)))

expr61_sub <- expr61_final[common_genes, ]
expr62_sub <- expr62_final[common_genes, ]
expr63_sub <- expr63_final[common_genes, ]

master_expr <- cbind(expr61_sub, expr62_sub, expr63_sub)

# Process and clean metadata for all three studies
meta61_f <- res61$metadata %>% 
  mutate(Diagnosis = clean_dx(`category:ch1`), Study = "GSE73461") %>%
  dplyr::select(title, Diagnosis, Study) # Explicitly use dplyr

meta62_f <- res62$metadata %>% 
  mutate(Diagnosis = clean_dx(`category:ch1`), Study = "GSE73462") %>%
  dplyr::select(title, Diagnosis, Study)

meta63_f <- res63$metadata %>% 
  mutate(Diagnosis = clean_dx(`category:ch1`), Study = "GSE73463") %>%
  dplyr::select(title, Diagnosis, Study)

master_metadata_combined <- rbind(meta61_f, meta62_f, meta63_f)

master_metadata_final <- master_metadata_combined %>%
  filter(title %in% colnames(master_expr)) %>%
  arrange(match(title, colnames(master_expr)))

design <- model.matrix(~Diagnosis, data = master_metadata_final)

master_expr_clean <- removeBatchEffect(
  master_expr, 
  batch = master_metadata_final$Study,
  design = design
)

set.seed(42)

umap_raw <- umap(t(master_expr))

plot_df_raw <- data.frame(
  UMAP1 = umap_raw$layout[,1],
  UMAP2 = umap_raw$layout[,2],
  Study = master_metadata_final$Study,
  Diagnosis = master_metadata_final$Diagnosis
)

ggplot(plot_df_raw, aes(x = UMAP1, y = UMAP2, color = Study, shape = Diagnosis)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP: Before Correction (Raw Merged)",
       subtitle = "Expect samples to cluster primarily by Study")


set.seed(42)
umap_clean <- umap(t(master_expr_clean))

plot_df_clean <- data.frame(
  UMAP1 = umap_clean$layout[,1],
  UMAP2 = umap_clean$layout[,2],
  Study = master_metadata_final$Study,
  Diagnosis = master_metadata_final$Diagnosis
)

# 3. Plot Cleaned
ggplot(plot_df_clean, aes(x = UMAP1, y = UMAP2, color = Study, shape = Diagnosis)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP: After Limma removeBatchEffect",
       subtitle = "Expect KD samples from all studies to overlap")

outlier_indices <- as.numeric(rownames(plot_df_clean[plot_df_clean$UMAP1 > 15, ]))
outlier_info <- master_metadata_final[outlier_indices, ]
print(outlier_info)

group <- factor(master_metadata_final$Diagnosis, levels = c("Control", "KD"))
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "KD_vs_Control")

fit <- lmFit(master_expr_clean, design)
fit <- eBayes(fit)

res_de <- topTable(fit, coef = "KD_vs_Control", number = Inf, sort.by = "P")

print(head(res_de, 10))

gene_key <- rbind(annov4[, c("ENTREZID", "SYMBOL")], 
                  annov3[, c("ENTREZID", "SYMBOL")]) %>%
  distinct(ENTREZID, .keep_all = TRUE)
res_de$ENTREZID <- rownames(res_de)
res_de_symbols <- merge(res_de, gene_key, by = "ENTREZID", all.x = TRUE)
res_de_symbols <- res_de_symbols[order(res_de_symbols$P.Value), ]
print(head(res_de_symbols[, c("SYMBOL", "logFC", "adj.P.Val", "ENTREZID")], 20))

entrez_ranks <- res_de$t
names(entrez_ranks) <- rownames(res_de) # These are your Entrez IDs
entrez_ranks <- na.omit(entrez_ranks)
entrez_ranks <- sort(entrez_ranks, decreasing = TRUE)

all_gene_sets <- msigdbr(species = "human", collection = "H")
pathways_entrez <- split(x = as.character(all_gene_sets$ncbi_gene), 
                         f = all_gene_sets$gs_name)
pathways_entrez <- lapply(pathways_entrez, function(x) x[!is.na(x) & x != ""])
names(entrez_ranks) <- as.character(names(entrez_ranks))
set.seed(42)
fgsea_entrez <- fgsea(pathways = pathways_entrez, 
                      stats = entrez_ranks,
                      minSize = 15,
                      maxSize = 500)

kd_indices <- which(master_metadata_final$Diagnosis == "KD")
kd_expr <- master_expr_clean[, kd_indices]
gene_vars <- apply(kd_expr, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
kd_expr_subset <- kd_expr[top_genes, ]
kd_expr_subset <- na.omit(kd_expr_subset)
rcc_results <- RecursiveConsensusClustering:::ccRun(
  d = kd_expr_subset, 
  maxK = 6, 
  repCount = 100, 
  clusterAlg = "km", 
  distance = "euclidean",
  pItem = 0.8,         # Standard: subsample 80% of items
  pFeature = 1,        # Use all features in subsamples
  verbose = TRUE
)
kd_matrix <- as.matrix(kd_expr_subset)
class(kd_matrix) <- "matrix"
kd_dist <- dist(t(kd_matrix), method = "euclidean")
rcc_plots <- ConsensusClusterPlus(
  d = kd_dist,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "km",
  distance = "euclidean",
  title = "KD_Recursive_Clustering",
  plot = "pdf",
  writeTable = TRUE
)

kd_subgroups <- rcc_plots[[3]][["consensusClass"]]

kd_meta <- master_metadata_final[kd_indices, ]
kd_meta$Subgroup <- as.factor(kd_subgroups)

table(kd_meta$Subgroup)

study_table <- table(kd_meta$Subgroup, kd_meta$Study)
print(study_table)
chisq.test(study_table)

meta_batch <- data.frame(study = master_metadata_final$Study[kd_indices])
harm_out <- harmony::HarmonyMatrix(
  data_mat = t(kd_matrix), 
  meta_data = meta_batch, 
  vars_use = 'study', 
  do_pca = TRUE
)

kd_dist_harmony <- dist(harm_out)
rcc_harmony <- ConsensusClusterPlus(
  d = kd_dist_harmony,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  clusterAlg = "km",
  distance = "euclidean",
  title = "KD_Harmony_Integrated",
  plot = "pdf"
)