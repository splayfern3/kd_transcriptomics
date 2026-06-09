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
library(sigclust)
library(pheatmap)
library(plotly)
library(patchwork)
library(ggforce)
library(scales)
library(GSVA)
library(GEOquery)
library(preprocessCore)
library(ggrepel)
library(patchwork)
source("r_code/illumina_loader.R")
source("r_code/gene_collapse.R")
source("r_code/meta_filter.R")
source("r_code/process_study.R")
source("r_code/clean_dx.R")

install.packages(c("devtools", "clv", "fields", "matrixStats", "data.table", "cluster", "clue", "circlize", "gdata"))


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")


devtools::install_github("MSCTR/RecursiveConsensusClustering")
library(RecursiveConsensusClustering)
                         
# TODO: add getGEO w/ conditionals so pipeline can be ran w/o assuming files already downloaded
# gse73463_meta <- getGEO("GSE73463", destdir = "transcriptome_data/uncompressed", getGPL = FALSE)




#load gse73461 to gse73463 along with the expression and pvalue files
gse73461 <- load_illumina_geotxt(
  path = "transcriptome_data/uncompressed/GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt"
)
expr73461_lin <- gse73461$expr_lin
pval73461 <- gse73461$pval
#Remove columns filled with NA
na_cols61 <- which(colSums(is.na(expr73461_lin)) == nrow(expr73461_lin))
if(length(na_cols61) > 0) expr73461_lin <- expr73461_lin[, -na_cols61]
expr73461_lin[expr73461_lin < 1] <- 1
expr73461_log2 <- log2(expr73461_lin)
expr73461_log2 <- normalize.quantiles(as.matrix(expr73461_log2))
rownames(expr73461_log2) <- rownames(expr73461_lin)
colnames(expr73461_log2) <- colnames(expr73461_lin)
rm(expr73461_lin); gc()
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

raw_81 <- read.delim("transcriptome_data/compressed/GSE63881_Non-normalized_Data.txt.gz",
                     header = TRUE, sep = "\t")
gse63881_meta_obj <- getGEO("GSE63881", destdir = "transcriptome_data/uncompressed", getGPL = FALSE)
gse63881_meta <- pData(gse63881_meta_obj[[1]])
intensity_cols81 <- grep("^X.*(?<!Detection\\.Pval)$", colnames(raw_81), perl=TRUE, value=TRUE)
pval_cols81 <- grep("Detection\\.Pval", colnames(raw_81), value=TRUE)
expr81_lin <- raw_81[, intensity_cols81]
rownames(expr81_lin) <- raw_81$ID_REF
pval81 <- raw_81[, pval_cols81]
rownames(pval81) <- raw_81$ID_REF
colnames(expr81_lin) <- gsub("^X", "", colnames(expr81_lin))
colnames(pval81) <- colnames(expr81_lin)
expr81_lin[expr81_lin < 1] <- 1
expr81_log2 <- log2(expr81_lin)
expr81_log2 <- normalize.quantiles(as.matrix(expr81_log2))
rownames(expr81_log2) <- rownames(expr81_lin)
colnames(expr81_log2) <- colnames(expr81_lin)
rm(expr81_lin); gc()

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


# For GSE73463
colnames(expr73463_log2) <- gsub("_Detection_Pval|_Detection Pval", "", colnames(expr73463_log2))
colnames(pval73463)      <- gsub("_Detection_Pval|_Detection Pval", "", colnames(pval73463))
meta63_full <- pData(gse73463_meta[[1]])
# Ensure rownames match title
rownames(meta63_full) <- meta63_full$title
# GSE63881 was already cleaned above
# filter out rows with high p values (i.e their intensity is comparable to the background)

# --- Process GSE73461 ---
res61 <- filter_by_metadata("GSE73461", expr73461_log2, pval73461, threshold = 0.5)
res61$metadata <- res61_meta_clean[intersect(colnames(res61$expr), rownames(res61_meta_clean)), ]

# --- Process GSE73462 ---
res62 <- filter_by_metadata("GSE73462", expr73462_log2, pval73462, threshold = 0.5)
res62$metadata <- meta62_full[intersect(colnames(res62$expr), rownames(meta62_full)), ]

# --- Process GSE73463 ---
res63 <- filter_by_metadata("GSE73463", expr73463_log2, pval73463, threshold = 0.5)
res63$metadata <- meta63_full[intersect(colnames(res63$expr), rownames(meta63_full)), ]

# --- Process GSE63881 ---
pval81 <- pval81[complete.cases(pval81), ] #Remove NAs
expr81_log2 <- expr81_log2[rownames(pval81), ]
mapping_table <- data.frame(
  GSM_ID = rownames(gse63881_meta),
  Scan_ID = as.character(gse63881_meta$description),
  stringsAsFactors = FALSE
)


kd_phenotypes <- c("Mild", "Moderate", "Severe")
target_phases <- c("Acute", "Convalescent")


meta81_kd_all <- gse63881_meta %>%
  filter(`phenotype:ch1` %in% kd_phenotypes) %>%
  filter(`phase:ch1` %in% target_phases)

# Align the Expression and P-value matrices
common_samples <- intersect(colnames(expr81_log2), rownames(meta81_kd_all))
expr81_final_subset <- expr81_log2[, common_samples]
pval81_final_subset <- pval81[, common_samples]

# Perform the p-value filtering
# Probes must be detected (p < 0.05) in more than 50% of these samples
keep_probes <- rowMeans(pval81_final_subset < 0.05, na.rm = TRUE) > 0.5
expr81_filtered <- expr81_final_subset[keep_probes, ]

#Create the result object
res81 <- list(
  expr = expr81_filtered,
  metadata = meta81_kd_all[common_samples, ]
)


gplv4 <- read.delim("transcriptome_data/platforms/v4/GPL10558.txt", comment.char = "#", stringsAsFactors = FALSE)
gplv3 <- read.delim("transcriptome_data/platforms/v3/GPL6947.txt", comment.char = "#", stringsAsFactors = FALSE)
#map probe ids to gene name and symbol 
annov4 <- gplv4[, c("Array_Address_Id", "ID", "Entrez_Gene_ID", "Symbol")]
annov3 <- gplv3[, c("Array_Address_Id", "Entrez_Gene_ID", "Symbol")]
#rename cols
colnames(annov4) <- c("PROBEID", "ID", "ENTREZID", "SYMBOL")
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
mapped81_for_collapse <- annov4 %>% 
  filter(ID %in% rownames(res81$expr)) %>% 
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

# GSE63881 was already filtered for KD samples because it uses different naming schemes than the other 3 studies, 
# thus processing must be done without the process_study function used above

mapped_81 <- mapped81_for_collapse[mapped81_for_collapse$ID %in% rownames(res81$expr), ]
expr81_ready <- res81$expr[as.character(mapped_81$ID), ]
expr81_collapsed <- WGCNA::collapseRows(
  datET = expr81_ready,
  rowGroup = mapped_81$ENTREZID,
  rowID = mapped_81$ID,
  method = "MaxMean"
)$datETcollapsed

expr81_final <- as.matrix(expr81_collapsed)

# Find genes present in all 4 datasets
common_all <- intersect(rownames(expr61_final), 
                        intersect(rownames(expr62_final), 
                                  intersect(rownames(expr63_final), rownames(expr81_final))))

message("Final gene count for merge: ", length(common_all))
expr61_sub <- expr61_final[common_all, ]
expr62_sub <- expr62_final[common_all, ]
expr63_sub <- expr63_final[common_all, ]
expr81_sub <- expr81_final[common_all, ]
master_expr <- cbind(expr61_sub, expr62_sub, expr63_sub, expr81_sub)

# Process and clean metadata for all 4 studies
# Process and clean metadata for all 4 studies with Gender and Acuity, Age, IVIG status, and Aneursym status
meta61_f <- res61$metadata %>% 
  mutate(
    # 1. Add prefix to title
    title = paste0("GSE73461_", title),
    Diagnosis = clean_dx(`category:ch1`), 
    Severity = Diagnosis,
    Study = "GSE73461",
    Gender = gsub("gender: ", "", `characteristics_ch1.2`),
    Acuity = "Not Reported",
    Age_Months = NA,
    IVIG_Status = "Not reported",
    Aneurysm_Status = "Not Reported"
  ) %>%
  dplyr::select(title, Diagnosis, Severity, Study, Gender, Acuity, Age_Months, IVIG_Status, Aneurysm_Status)
# Also prefix the expression matrix columns
colnames(expr61_sub) <- paste0("GSE73461_", colnames(expr61_sub))

meta62_f <- res62$metadata %>% 
  mutate(
    title = paste0("GSE73462_", title),
    Diagnosis = clean_dx(`category:ch1`), 
    Severity = Diagnosis,
    Study = "GSE73462",
    Gender = gsub("gender: ", "", `characteristics_ch1.2`),
    Acuity = "Not Reported",
    Age_Months = NA,
    IVIG_Status = "Not Reported",
    Aneurysm_Status = "Not Reported"
  ) %>%
  filter(Diagnosis != "Uncertain") %>%
  dplyr::select(title, Diagnosis, Severity, Study, Gender, Acuity, Age_Months, IVIG_Status, Aneurysm_Status)

colnames(expr62_sub) <- paste0("GSE73462_", colnames(expr62_sub))

meta63_f <- res63$metadata %>% 
  mutate(
    title = paste0("GSE73463_", title),
    Diagnosis = clean_dx(`category:ch1`), 
    Severity = Diagnosis,
    Study = "GSE73463",
    Gender = gsub("gender: ", "", `characteristics_ch1.2`),
    Acuity = case_when(
      grepl("acute", `characteristics_ch1.1`, ignore.case = TRUE) ~ "Acute",
      grepl("conv", `characteristics_ch1.1`, ignore.case = TRUE) ~ "Convalescent",
      TRUE ~ "Unknown"
    ),
    Age_Months = NA,
    IVIG_Status = "Not Reported",
    Aneurysm_Status = "Not Reported"
  ) %>%
  dplyr::select(title, Diagnosis, Severity, Study, Gender, Acuity, Age_Months, IVIG_Status, Aneurysm_Status)

colnames(expr63_sub) <- paste0("GSE73463_", colnames(expr63_sub))

# Process GSE63881 metadata to match your master format
meta81_f <- res81$metadata %>%
  # Exclude OFI as requested
  filter(`phenotype:ch1` != "OFI") %>%
  mutate(
    title = paste0("GSE63881_", rownames(.)),
    Diagnosis = "KD", #all samples in 81 are KD since we filtered out OFI (other febral illness)
    Severity = `phenotype:ch1`,
    Study = "GSE63881",
    Gender = ifelse(`Sex:ch1` == 1, "male", "female"),
    Acuity = `phase:ch1`,
    Age_Months = as.numeric(`age_month:ch1`),
    IVIG_Status = `ivig:ch1`,
    Aneurysm_Status = `aneurysm:ch1`
  ) %>%
  dplyr::select(title, Diagnosis, Severity, Study, Gender, Acuity, Age_Months, IVIG_Status, Aneurysm_Status)

# Update column names in the expression matrix to match these titles
colnames(expr81_sub) <- paste0("GSE63881_", colnames(expr81_sub))


all_meta_list <- list(meta61_f, meta62_f, meta63_f, meta81_f)
master_metadata_combined <- bind_rows(all_meta_list)
master_expr <- cbind(expr61_sub, expr62_sub, expr63_sub, expr81_sub)
common_samples <- intersect(master_metadata_combined$title, colnames(master_expr))


master_metadata_final <- master_metadata_combined %>%
  filter(title %in% common_samples) %>%
  distinct(title, .keep_all = TRUE) %>%
  arrange(match(title, colnames(master_expr)))

rownames(master_metadata_final) <- master_metadata_final$title
master_expr_final <- master_expr[, master_metadata_final$title]


cb_palette <- c("GSE73461" = "#E69F00", # Orange
                "GSE73462" = "#56B4E9", # Sky Blue
                "GSE73463" = "#009E73", # Bluish Green
                "GSE63881" = "#CC79A7") # Reddish Purple




#####################
# ---  Analysis  ---#
#####################
# combat to remove batch effects between studies
mod <- model.matrix(~as.factor(Diagnosis), data = master_metadata_final)
master_expr_combat <- ComBat(dat = master_expr_final,
                             batch = master_metadata_final$Study, 
                             mod = mod, par.prior = TRUE)

kd_indices <- which(master_metadata_final$Diagnosis == "KD")
kd_expr <- master_expr_combat[, kd_indices]
kd_meta <- master_metadata_final[kd_indices, ]
kd_meta <- kd_meta %>%
  mutate(
    # Severity Mapping: Use specific status for 81, default to "KD" for others
    Aneurysm_Status = case_when(
      tolower(Aneurysm_Status) == "aneurysm" ~ "Aneurysm",
      tolower(Aneurysm_Status) == "dilated"  ~ "Dilated",
      tolower(Aneurysm_Status) == "normal"   ~ "Normal",
      TRUE ~ "KD"  # Sets Studies 61, 62, 63 as the baseline
    ),
    
    # IVIG Mapping: Use specific status for 81, default to "KD" for others
    IVIG_Status = case_when(
      tolower(IVIG_Status) == "responsive" ~ "Responsive",
      tolower(IVIG_Status) == "resistant"  ~ "Resistant",
      TRUE ~ "KD"  # Unified background label
    )
  )
gene_mads <- apply(kd_expr, 1, mad)
top_5k_genes <- names(sort(gene_mads, decreasing = TRUE))[1:5000]
kd_scaled <- t(scale(t(kd_expr[top_5k_genes, ])))
d_dist_scaled <- as.dist(1 - cor(kd_scaled, method = "pearson"))

# recursive clustering up to k = 10
set.seed(42)
results_k10 <- ConsensusClusterPlus(
  d = d_dist_scaled, 
  maxK = 10, #max of 10 clusters
  reps = 1000, # repeat 1000 times
  pItem = 0.8, # 80 percent of kd patients choosen
  clusterAlg = "hc",  
  innerLinkage = "ward.D2", 
  title = "kd_clusters", 
  plot = "png" 
)

# Assign the stable K=4 Subgroups
kd_meta$Subgroup <- as.factor(results_k10[[4]][["consensusClass"]])

# --- 3. SAVE THE FRESH PLOTS TO YOUR MAIN FOLDER ---
# Function to generate PCA dataframes
get_pca_df <- function(expr_matrix, meta) {
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Study = meta$Study)
  return(df)
}

# 2. Get PCA from Raw Merged Data (Before)
df_raw <- get_pca_df(master_expr_final, master_metadata_final)

# 3. Get PCA from ComBat Corrected Data (After)
df_combat <- get_pca_df(master_expr_combat, master_metadata_final)

# 4. Visualization
p_raw <- ggplot(df_raw, aes(x=PC1, y=PC2, color=Study)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw() + labs(title = "Raw Merged Data (No Correction)") +
  scale_color_manual(values = cb_palette)

p_combat <- ggplot(df_combat, aes(x=PC1, y=PC2, color=Study)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw() + labs(title = "Post-ComBat Correction") +
  scale_color_manual(values = cb_palette)

(p_raw | p_combat) + plot_layout(guides = 'collect')


# PCA Plot colored by subgroup, shape indicating acuity
pca_res <- prcomp(t(kd_scaled), scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], 
  Acuity = kd_meta$Acuity, Subgroup = kd_meta$Subgroup
)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subgroup, shape = Acuity)) +
  geom_point(alpha = 0.8, size = 5, stroke = 1.2) +
  theme_bw() +
  scale_color_manual(values = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#000B37")) +
  labs(title = "Final Subgroup Discovery (K=4)", subtitle = "Validated by PAC and Silhouette Analysis")
print(p_pca)

# PCA dataframe with clinical data
kd_shapes <- c("1" = 16, "2" = 17, "3" = 15, "4" = 3)
pca_df_clinical <- data.frame(
  PC1 = pca_res$x[,1], 
  PC2 = pca_res$x[,2], 
  Study = kd_meta$Study,
  Subgroup = as.factor(kd_meta$Subgroup),
  Severity = as.factor(kd_meta$Severity),     
  Aneurysm = as.factor(kd_meta$Aneurysm_Status), 
  IVIG = as.factor(kd_meta$IVIG_Status)         
)
# Plot A: Phenotype Classification
p1 <- ggplot(pca_df_clinical, aes(x = PC1, y = PC2, color = Severity, shape = Subgroup)) +
  geom_point(alpha = 0.8, size = 4) +
  scale_shape_manual(values = kd_shapes) +
  scale_color_manual(values = c("Severe"="#D55E00", "Moderate"="#E69F00", "Mild"="#009E73", "KD"="#999999")) +
  theme_bw() +
  labs(title = "Phenotype Classification", color = "Severity")

# Plot B: Aneurysm Status (Bottom Left)
p2 <- ggplot(pca_df_clinical, aes(x = PC1, y = PC2, color = Aneurysm, shape = Subgroup)) +
  geom_point(alpha = 0.8, size = 3.5, stroke = 1) +
  scale_shape_manual(values = kd_shapes) + # Set to plus sign for 4
  scale_color_manual(values = c("Aneurysm"="#D55E00", "Dilated"="#E69F00", "Normal"="#56B4E9", "KD"="#999999")) +
  theme_bw() +
  labs(title = "Clinical Severity: Aneurysm Status", color = "Aneurysm Status")

# Plot C: IVIG Response (Bottom Right)
p3 <- ggplot(pca_df_clinical, aes(x = PC1, y = PC2, color = IVIG, shape = Subgroup)) +
  geom_point(alpha = 0.8, size = 3.5, stroke = 1) +
  scale_shape_manual(values = kd_shapes) + # Set to plus sign for 4
  scale_color_manual(values = c("Resistant"="#CC79A7", "Responsive"="#E69F00", "KD"="#999999")) +
  theme_bw() +
  labs(title = "Clinical Severity: IVIG Response", color = "IVIG Status")

# --- STEP 4: Merge with Patchwork ---
combined_clinical_pca <- (p1) / (p2 | p3) + 
  plot_annotation(
    title = 'Integrative Clinical & Molecular PCA',
    
  )

print(combined_clinical_pca)

# --- STEP 1: Filter for Known Aneurysm Status ---
# We exclude "KD" (Not Reported) to match the other "Known Outcomes" plots
aneurysm_known <- kd_meta[kd_meta$Aneurysm_Status %in% c("Aneurysm", "Dilated", "Normal"), ]

# --- STEP 2: Statistical Test ---
# Using Fisher's with simulation to handle the "Zero" in Cluster 2
aneurysm_table <- table(aneurysm_known$Subgroup, aneurysm_known$Aneurysm_Status)
fisher_aneurysm <- fisher.test(aneurysm_table, simulate.p.value = TRUE, B = 10000)

# --- STEP 3: Data Prep for Proportional Plotting ---
df_aneurysm <- as.data.frame(aneurysm_table)
colnames(df_aneurysm) <- c("Subgroup", "Status", "Count")

df_aneurysm <- df_aneurysm %>%
  group_by(Subgroup) %>%
  mutate(Percent = Count / sum(Count))

# --- STEP 4: Generate the Plot ---
ggplot(df_aneurysm, aes(x = Subgroup, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white") +
  # Add white percentage labels
  geom_text(aes(label = ifelse(Percent > 0.05, scales::percent(Percent, accuracy = 1), "")), 
            position = position_fill(vjust = 0.5), color = "white", fontface = "bold") +
  # Using the Aneurysm-specific palette from your earlier PCA
  scale_fill_manual(values = c("Aneurysm" = "#D55E00", "Dilated" = "#E69F00", "Normal" = "#56B4E9")) +
  theme_classic(base_size = 14) +
  labs(title = "Coronary Risk by Subgroup",
       subtitle = paste0("Fisher's Exact p = ", format.pval(fisher_aneurysm$p.value)),
       x = "Molecular Cluster", y = "Proportion", fill = "Status")


# --- STEP 1: Filter for Known IVIG Status ---
ivig_known <- kd_meta[kd_meta$IVIG_Status %in% c("Resistant", "Responsive"), ]

# --- STEP 2: Statistical Test (Fisher's Exact) ---
ivig_table <- table(ivig_known$Subgroup, ivig_known$IVIG_Status)
# Using simulate.p.value because of sparse cells in some clusters
fisher_ivig <- fisher.test(ivig_table, simulate.p.value = TRUE, B = 10000)

# --- STEP 3: Plotting ---
df_ivig <- as.data.frame(ivig_table)
colnames(df_ivig) <- c("Subgroup", "Status", "Count")

df_ivig <- df_ivig %>%
  group_by(Subgroup) %>%
  mutate(Percent = Count / sum(Count))

ggplot(df_ivig, aes(x = Subgroup, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white", size = 0.2) +
  geom_text(aes(label = ifelse(Percent > 0.05, scales::percent(Percent, accuracy = 1), "")), 
            position = position_fill(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Resistant" = "#D55E00", "Responsive" = "#0072B2")) +
  theme_classic(base_size = 14) +
  labs(title = "IVIG Response by Subgroup",
       subtitle = paste0("Fisher's Exact p = ", format.pval(fisher_ivig$p.value)),
       x = "Molecular Cluster", y = "Proportion", fill = "Response")

# --- STEP 1: Filter for Known Severity ---
sev_known <- kd_meta[kd_meta$Severity %in% c("Mild", "Moderate", "Severe"), ]

# --- STEP 2: Statistical Test ---
sev_table <- table(sev_known$Subgroup, sev_known$Severity)
fisher_sev <- fisher.test(sev_table, simulate.p.value = TRUE, B = 10000)

# --- STEP 3: Plotting ---
df_sev <- as.data.frame(sev_table)
colnames(df_sev) <- c("Subgroup", "Level", "Count")

df_sev <- df_sev %>%
  group_by(Subgroup) %>%
  mutate(Percent = Count / sum(Count))

ggplot(df_sev, aes(x = Subgroup, y = Count, fill = Level)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white") +
  geom_text(aes(label = ifelse(Percent > 0.05, scales::percent(Percent, accuracy = 1), "")), 
            position = position_fill(vjust = 0.5), color = "white", fontface = "bold") +
  # Nature-style diverging palette for severity
  scale_fill_manual(values = c("Mild" = "#009E73", "Moderate" = "#E69F00", "Severe" = "#D55E00")) +
  theme_classic(base_size = 14) +
  labs(title = "Clinical Severity Distribution",
       subtitle = paste0("Fisher's Exact p = ", format.pval(fisher_sev$p.value)),
       x = "Molecular Cluster", y = "Proportion", fill = "Severity")




















# --- STEP 1: Define the Model (The "Discovery" Engine) ---
group_list <- factor(kd_meta$Subgroup)
design <- model.matrix(~0 + group_list)
colnames(design) <- c("C1", "C2", "C3", "C4")

fit <- lmFit(kd_expr, design)

# Create the One-vs-Rest contrasts
contrast_matrix <- makeContrasts(
  C1_vs_All = C1 - (C2 + C3 + C4)/3,
  C2_vs_All = C2 - (C1 + C3 + C4)/3,
  C3_vs_All = C3 - (C1 + C2 + C4)/3,
  C4_vs_All = C4 - (C1 + C2 + C3)/3,
  levels = design
)

fit_clus2 <- contrasts.fit(fit, contrast_matrix)
fit_clus2 <- eBayes(fit_clus2)

# --- STEP 2: Use fit_clus2 to get your Heatmap Genes ---
# This pulls the 50 genes with the highest F-statistic (most variation across groups)
top_table_all <- topTable(fit_clus2, number = 50, sort.by = "F")
top_50_safe <- rownames(top_table_all)

# Now your heatmap code below will work!
matrix_ids <- rownames(kd_scaled[top_50_safe, ])
mapped_symbols <- mapIds(org.Hs.eg.db, keys = matrix_ids, column = "SYMBOL", 
                         keytype = "ENTREZID", multiVals = "first")
final_labels <- ifelse(!is.na(mapped_symbols), as.character(mapped_symbols), matrix_ids)

# --- STEP 3: Define a Hierarchical Column Order ---
# This sorts samples by Subgroup (1-4), then by IVIG_Status, then by CAA_Status.
# Replace "IVIG_Status" and "Acuity" with the exact column names in your kd_meta.
final_order <- order(
  kd_meta$Subgroup, 
  kd_meta$IVIG_Status, 
  kd_meta$Aneurysm_Status # Sorting by CAA/Aneurysm status within Subgroups
)

# --- STEP 4: Define Annotation Colors for Clinical Outcomes ---
# This ensures the new clinical rows have professional, distinct colors.
anno_colors <- list(
  Subgroup = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#000B37"),
  Acuity = c("Acute" = "#FF4500", "Convalescent" = "#1E90FF", "Not Reported" = "#999999"),
  IVIG_Status = c("Responsive" = "#0072B2", "Resistant" = "#FF8300", "KD" = "#999999"),
  Aneurysm_Status = c("Normal" = "#56B4E9", "Dilated" = "#E69F00", "Aneurysm" = "#CDFFCD", "KD" = "#999999")
)

# --- STEP 5: Plot the Heatmap ---
pheatmap(
  kd_scaled[top_50_safe, final_order], 
  labels_row = final_labels, 
  show_colnames = FALSE, 
  cluster_cols = FALSE, # Required to maintain our custom sort order
  cluster_rows = TRUE,  # Keeps similar genes grouped together
  
  # Adding IVIG and Aneurysm status to the top annotations
  annotation_col = kd_meta[final_order, c("Acuity", "IVIG_Status", "Aneurysm_Status", "Subgroup")],
  annotation_colors = anno_colors,
  
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Top 50 Gene Markers Organized by Subgroup and Clinical Outcome"
)
# Force plots to appear in the RStudio Plots pane
p1 <- ggplot(val_df, aes(x=K, y=Silhouette)) +
  geom_line(color="steelblue", linewidth=1) + 
  geom_point(size=3) +
  theme_bw() + 
  labs(title="Silhouette Score", subtitle="Higher = Better Separation", y="Mean Width")

p2 <- ggplot(val_df, aes(x=K, y=PAC)) +
  geom_line(color="firebrick", linewidth=1) + 
  geom_point(size=3) +
  theme_bw() + 
  labs(title="PAC Score", subtitle="Lower = Most Stable Clusters", y="Ambiguity Index")

# This triggers the display
p1 + p2


# Function for standardized, professional volcano plots
make_poster_volcano <- function(fit_obj, coef_name, cluster_num) {
  df <- topTable(fit_obj, coef = coef_name, number = Inf)
  df$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(df), column = "SYMBOL", keytype = "ENTREZID")
  
  # STATISTICAL LIMITS (Adjust these to change what gets colored)
  p_thresh <- 0.05    # FDR threshold (adj.P.Val)
  fc_thresh <- 1.0    # X-axis threshold (2-fold change)

  # Define coloring logic
  df$diffexpressed <- "Not Significant"
  df$diffexpressed[df$logFC > fc_thresh & df$adj.P.Val < p_thresh] <- "Upregulated"
  df$diffexpressed[df$logFC < -fc_thresh & df$adj.P.Val < p_thresh] <- "Downregulated"
  
  # Get top 20 genes to label
  top_labels <- rbind(
    head(df[df$diffexpressed == "Upregulated", ], 10),
    head(df[df$diffexpressed == "Downregulated", ], 10)
  )
  
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.3, size = 1.2) +
    # STANDARDIZED COLORS: Red for Up, Blue for Down, Grey for No
    scale_color_manual(values = c("Downregulated" = "dodgerblue3", 
                                  "Not Significant" = "grey80", 
                                  "Upregulated" = "firebrick3")) +
    geom_text_repel(data = top_labels, aes(label = Symbol), 
                    size = 3.5, fontface = "bold", 
                    box.padding = 0.5, segment.color = 'grey50',
                    max.overlaps = Inf) +
    # Removed subtitle, kept title only
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 16, face = "bold")) +
    labs(title = paste("Cluster", cluster_num), 
         x = "Log2 Fold Change", 
         y = "-log10 P-value")
}

# Generate the 4 plots with consistent colors
v1 <- make_poster_volcano(fit_clus2, "C1_vs_All", "1")
v2 <- make_poster_volcano(fit_clus2, "C2_vs_All", "2")
v3 <- make_poster_volcano(fit_clus2, "C3_vs_All", "3")
v4 <- make_poster_volcano(fit_clus2, "C4_vs_All", "4")

# Combine into a 2x2 grid
(v1 | v2) / (v3 | v4)

# 1. Prepare Gene Sets (Updated for msigdbr 10.0.0+)
m_df <- msigdbr(species = "Homo sapiens", collection = "H") # 'H' is for Hallmark
pathways_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
gsea_results <- list()

# 2. Run the Analysis Loop
for(i in 1:4){
  coef_name <- paste0("C", i, "_vs_All")
  tab <- topTable(fit_clus2, coef = coef_name, number = Inf)
  
  # Mapping Symbols
  tab$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(tab), column = "SYMBOL", keytype = "ENTREZID")
  
  # Create Ranks
  ranks <- tab$t
  names(ranks) <- tab$Symbol
  
  # CLEANING: Remove NAs, keep highest |t| for duplicated symbols
  ranks <- ranks[!is.na(names(ranks))]
  ranks_df <- data.frame(symbol = names(ranks), t = ranks)
  ranks_df <- ranks_df[order(-abs(ranks_df$t)), ]
  ranks_df <- ranks_df[!duplicated(ranks_df$symbol), ]
  ranks <- setNames(ranks_df$t, ranks_df$symbol)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Run GSEA
  gsea_results[[paste0("C", i)]] <- fgsea(pathways = pathways_list, 
                                          stats = ranks, 
                                          minSize=15, 
                                          maxSize=500)
}

# 3. Create the Summary Table
summary_df <- data.frame(
  Cluster = 1:4,
  Top_Pathway = sapply(gsea_results, function(x) x[order(padj)][1, pathway]),
  NES = sapply(gsea_results, function(x) round(x[order(padj)][1, NES], 2)),
  P_Adj = sapply(gsea_results, function(x) format.pval(x[order(padj)][1, padj], digits=3))
)

# 1. Collect top 2 pathways for each cluster to show more data
plot_data <- bind_rows(
  lapply(names(gsea_results), function(name) {
    gsea_results[[name]] %>%
      arrange(padj) %>%
      head(2) %>%
      mutate(Cluster = name)
  })
)
# 2. Clean up names for the poster
plot_data$pathway <- gsub("HALLMARK_", "", plot_data$pathway)
plot_data$pathway <- gsub("_", " ", plot_data$pathway)

# 3. Create the Directional Plot
ggplot(plot_data, aes(x = NES, y = reorder(pathway, NES), fill = NES > 0)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "dodgerblue3"), 
                    labels = c("Downregulated", "Upregulated")) +
  facet_wrap(~Cluster, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold")) +
  labs(title = "Top Biological Signatures by Cluster",
       x = "Normalized Enrichment Score (NES)",
       y = NULL,
       fill = "Direction")



# --- STEP 1: DATA MERGING & PREFIX CLEANING ---

hc_samples <- master_metadata_final[master_metadata_final$Diagnosis == "Control", ]
hc_samples$Subgroup <- "control"

kd_to_join <- kd_meta[, c("title", "Subgroup")]
hc_to_join <- hc_samples[, c("title", "Subgroup")]
combined_meta <- rbind(kd_to_join, hc_to_join)

# Strip GSE prefix to match matrix column names
combined_meta$clean_title <- sub("^GSE[0-9]+_", "", combined_meta$title)

common_samples <- intersect(combined_meta$clean_title, colnames(master_expr_combat))
message("Successfully Synced Samples: ", length(common_samples))

combined_meta_final <- combined_meta[combined_meta$clean_title %in% common_samples, ]
combined_expr_final <- master_expr_combat[, combined_meta_final$clean_title]

colnames(combined_expr_final) <- combined_meta_final$title
combined_meta_final$Group <- factor(paste0("Cluster_", combined_meta_final$Subgroup))


# --- STEP 2: STATISTICAL MODELING (LIMMA) ---

design_hc <- model.matrix(~0 + Group, data = combined_meta_final)
colnames(design_hc) <- make.names(colnames(design_hc))

fit_hc <- lmFit(combined_expr_final, design_hc)

contrast_hc <- makeContrasts(
  C1_vs_Ctrl = GroupCluster_1 - GroupCluster_control,
  C2_vs_Ctrl = GroupCluster_2 - GroupCluster_control,
  C3_vs_Ctrl = GroupCluster_3 - GroupCluster_control,
  C4_vs_Ctrl = GroupCluster_4 - GroupCluster_control,
  levels = design_hc
)

fit_hc2 <- contrasts.fit(fit_hc, contrast_hc)
fit_hc2 <- eBayes(fit_hc2)


# --- STEP 3: PATHWAY ENRICHMENT (GSEA) ---

pathways_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
gsea_hc_list <- list()

for(i in 1:4){
  coef_name <- paste0("C", i, "_vs_Ctrl")
  tab <- topTable(fit_hc2, coef = coef_name, number = Inf)
  
  # Map Entrez IDs to Gene Symbols
  tab$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(tab), 
                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  
  # --- THE FIX: Remove NAs from Symbols BEFORE naming ranks ---
  tab_clean <- tab[!is.na(tab$Symbol), ]
  
  # Rank genes by t-statistic
  ranks <- tab_clean$t
  names(ranks) <- tab_clean$Symbol
  
  # Remove duplicates (keep highest |t|) and sort
  ranks_df2 <- data.frame(symbol = names(ranks), t = ranks)
  ranks_df2 <- ranks_df2[order(-abs(ranks_df2$t)), ]
  ranks_df2 <- ranks_df2[!duplicated(ranks_df2$symbol), ]
  ranks <- sort(setNames(ranks_df2$t, ranks_df2$symbol), decreasing = TRUE)
  
  # Run GSEA
  gsea_hc_list[[paste0("Cluster ", i)]] <- fgsea(pathways = pathways_list, 
                                                 stats = ranks, 
                                                 minSize=15, 
                                                 maxSize=500)
}


# 1. Identify the Top 10 Up and Top 10 Down for each cluster
top_pathways_list <- lapply(names(gsea_results), function(name) {
  res <- gsea_results[[name]] %>% filter(padj < 0.05)
  
  up <- res %>% slice_max(order_by = NES, n = 10)
  down <- res %>% slice_min(order_by = NES, n = 10)
  
  return(unique(c(up$pathway, down$pathway)))
})

# 2. Get the unique set of all these pathways (The Union)
all_top_pathways <- unique(unlist(top_pathways_list))

# 3. Filter the full results for just these pathways
dot_plot_df <- bind_rows(lapply(names(gsea_results), function(name) {
  gsea_results[[name]] %>% 
    filter(pathway %in% all_top_pathways) %>%
    mutate(Cluster = name)
}))

# 4. Clean Pathway Names
dot_plot_df$pathway_clean <- gsub("_", " ", gsub("HALLMARK_", "", dot_plot_df$pathway))

# 5. Create the Dot Plot
ggplot(dot_plot_df, aes(x = Cluster, y = reorder(pathway_clean, NES))) +
  # Use size for magnitude and color for direction/significance
  geom_point(aes(size = abs(NES), color = NES)) +
  
  # Professional Red/Blue Divergent Scale
  scale_color_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0) +
  
  # Theme and Labels
  theme_bw() +
  labs(title = "Comparative Pathway Signatures: Top 10 Up/Down per Cluster",
       subtitle = "Dot size = Magnitude of Enrichment | Color = Direction (NES)",
       x = "Molecular Subgroup",
       y = "Hallmark Pathway",
       size = "abs(NES)",
       color = "NES Score") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "grey95"),
    legend.position = "right"
  )
# --- PCA: PC1 vs PC2 with Control Diamonds ---
pca_full <- prcomp(t(master_expr_combat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_full$x[,1], 
  PC2 = pca_full$x[,2], 
  Group = combined_meta_final$Group[match(colnames(master_expr_combat), combined_meta_final$title)]
)

p_pca_final <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(values = c("Cluster_1" = 16, "Cluster_2" = 16, 
                                "Cluster_3" = 16, "Cluster_4" = 16, 
                                "Cluster_control" = 18)) + 
  scale_color_manual(values = c("Cluster_1" = "#E69F00", "Cluster_2" = "#56B4E9", 
                                "Cluster_3" = "#009E73", "Cluster_4" = "#000B37", 
                                "Cluster_control" = "#D55E00")) +
  theme_minimal() +
  labs(title = "PCA: KD Subgroups & Healthy Controls")

print(p_pca_final)

# --- Standalone Violin Grid ---


p_violins_sina <- ggplot(le_plot_df, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(alpha = 0.2, scale = "width", color = "grey30", trim = TRUE) + 
  geom_sina(aes(color = Group), size = 1.2, alpha = 0.6, maxwidth = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black", alpha = 0.5) +
  # Fixed: changed ncol to 4 and adjusted text size to prevent title clipping
  facet_wrap(~Pathway, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Cluster_1" = "#E69F00", "Cluster_2" = "#56B4E9", 
                               "Cluster_3" = "#009E73", "Cluster_4" = "#000B37", 
                               "Cluster_control" = "#D55E00")) +
  scale_color_manual(values = c("Cluster_1" = "#E69F00", "Cluster_2" = "#56B4E9", 
                                "Cluster_3" = "#009E73", "Cluster_4" = "#000B37", 
                                "Cluster_control" = "#D55E00")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey95"),
    # Fixed: decreased size slightly to ensure "REACTIVE OXYGEN SPECIES" fits
    strip.text = element_text(face = "bold", size = 8), 
    legend.position = "right"
  ) +
  labs(title = "Leading Edge Pathway Scores with Individual Sample Distribution",
       y = "Z-score Normalized Expression", x = NULL)

print(p_violins_sina)




# 1. Setup Groups
group_list <- factor(kd_meta$Subgroup)
design <- model.matrix(~0 + group_list)
colnames(design) <- c("C1", "C2", "C3", "C4")

fit <- lmFit(kd_expr, design)

# 2. One-vs-Rest Contrasts
# This finds genes unique to EACH subgroup
contrast_matrix <- makeContrasts(
  C1_unique = C1 - (C2 + C3 + C4)/3,
  C2_unique = C2 - (C1 + C3 + C4)/3,
  C3_unique = C3 - (C1 + C2 + C4)/3,
  C4_unique = C4 - (C1 + C2 + C3)/3,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_markers <- lapply(1:4, function(i) {
  # Using number = Inf ensures we don't miss any genes
  # Removing sort.by allows us to handle the sorting in our master pipe
  topTable(fit2, coef = i, number = Inf) 
})
# Important: Name the list so the .id = "cluster" step works
names(all_markers) <- c("1", "2", "3", "4")

target_genes <- all_markers %>%
  lapply(function(x) rownames_to_column(as.data.frame(x), var = "gene")) %>%
  bind_rows(.id = "cluster") %>%
  filter(abs(logFC) > 1) %>%
  group_by(cluster) %>%
  slice_min(order_by = P.Value, n = 10, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()
unique_marker_ids <- target_genes
# 1. Pull these 40 unique genes from the FULL KD expression matrix
# This avoids the "out of bounds" error because we go back to the source
raw_marker_mat <- kd_expr[unique_marker_ids, ]

# 2. Scale them (Z-score) so we can see the differences in a heatmap
# Math: (Value - Mean) / SD
plot_mat <- t(scale(t(raw_marker_mat)))[, final_order]

# 3. Map the Entrez IDs (like 146723) to Gene Symbols
mapped_symbols <- mapIds(org.Hs.eg.db, 
                         keys = rownames(plot_mat), 
                         column = "SYMBOL", 
                         keytype = "ENTREZID", 
                         multiVals = "first")

# 4. Create row labels: Use Symbol if found, otherwise keep the ID
final_labels <- ifelse(!is.na(mapped_symbols), as.character(mapped_symbols), rownames(plot_mat))

# 5. Generate the Heatmap
pheatmap(plot_mat, 
         labels_row = final_labels, 
         show_colnames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, # Keep them grouped by cluster (10 genes each)
         annotation_col = kd_meta[final_order, c("Acuity", "IVIG_Status", "Aneurysm_Status", "Subgroup")],
         annotation_colors = anno_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 10 Markers per Cluster")

# 1. Create a contingency table (Counts of Acuity per Cluster)
clinical_table <- table(kd_meta$Subgroup, kd_meta$Acuity)
print("Contingency Table: Subgroup vs Acuity")
print(clinical_table)

# 2. Run the Chi-Squared Test
# Math: Compares 'Observed' counts to 'Expected' counts if groups were random
chi_test <- chisq.test(clinical_table)
print(chi_test)

# 3. Calculate 'Residuals' to see WHICH group is driving the significance
# Positive residuals = more samples than expected (Enrichment)
# Negative residuals = fewer samples than expected (Depletion)
print("Pearson Residuals (Strength of association):")
print(chi_test$residuals)

# 1. Fix the grouping function and calculate percentages
df_plot <- as.data.frame(clinical_table)
colnames(df_plot) <- c("Subgroup", "Acuity", "Count")

df_plot <- df_plot %>%
  group_by(Subgroup) %>%  # Fixed function name here
  mutate(Percent = Count / sum(Count))

# 2. Plot with "Nature Portfolio" style aesthetics
ggplot(df_plot, aes(x = Subgroup, y = Count, fill = Acuity)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "white", size = 0.2) +
  # Add bold white labels only for groups larger than 5%
  geom_text(aes(label = ifelse(Percent > 0.05, percent(Percent, accuracy = 1), "")), 
            position = position_fill(vjust = 0.5), 
            color = "white", fontface = "bold", size = 3.5) +
  # High-end color palette (Muted, professional tones)
  scale_fill_manual(values = c("Acute" = "#D55E00",      # Vermillion
                               "Convalescent" = "#0072B2", # Blue
                               "Not Reported" = "#999999")) + # Grey
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(size = 0.8),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Clinical Enrichment by Molecular Subgroup",
    subtitle = paste0("Chi-squared p = ", format.pval(chi_test$p.value)),
    x = "Molecular Cluster",
    y = "Proportion of Samples",
    fill = "Clinical State"
  )

gender_table <- table(kd_meta$Subgroup, kd_meta$Gender)
print("Contingency Table: Subgroup vs Gender")
print(gender_table)

library(GSVA)
library(msigdbr)

# 1. Download the "Hallmark" gene sets
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- split(x = h_df$entrez_gene, f = h_df$gs_name)

# 2. Setup the Parameter Object (Fixes the 'inherited method' error)
# Math: kcdf="Gaussian" is for log-normalized data; use "Poisson" for raw counts.
gsvapar <- gsvaParam(exprData = as.matrix(kd_expr), 
                     geneSets = h_list, 
                     kcdf = "Gaussian")

# 3. Run GSVA
gsva_results <- gsva(gsvapar)

# 1. Define extreme color breaks
# This makes any score > 0.4 deep red and < -0.4 deep blue
my_breaks <- seq(-0.4, 0.4, length.out = 101)

# 2. Re-run the plot with updated scaling
pheatmap(gsva_results[, final_order],
         show_colnames = FALSE,
         cluster_cols = FALSE,
         annotation_col = kd_meta[final_order, c("Acuity", "IVIG_Status", "Aneurysm_Status", "Subgroup")],
         annotation_colors = anno_colors,
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         breaks = my_breaks, # <--- Add this line for extreme scaling
         main = "Figure 3: Hallmark Pathway Variation (High Contrast)")

# Helper function to process and map symbols/descriptions
process_marker_list <- function(marker_list) {
  lapply(names(marker_list), function(cid) {
    df <- as.data.frame(marker_list[[cid]])
    df$Entrez_ID <- rownames(df)
    df$Cluster <- cid
    return(df)
  }) %>%
    bind_rows() %>%
    mutate(
      Symbol = mapIds(org.Hs.eg.db, keys=Entrez_ID, column="SYMBOL", keytype="ENTREZID", multiVals="first"),
      Description = mapIds(org.Hs.eg.db, keys=Entrez_ID, column="GENENAME", keytype="ENTREZID", multiVals="first")
    ) %>%
    select(Cluster, Symbol, Entrez_ID, Description, logFC, P.Value, adj.P.Val) %>%
    arrange(Cluster, P.Value)
}

# --- GENERATE TABLE 1: Cluster vs. Other KD (Subgroup Discovery) ---
table_subgroups <- process_marker_list(all_markers) # Assumes all_markers came from fit_clus2
write.csv(table_subgroups, "Table1_Cluster_vs_Cluster_Full.csv", row.names = FALSE)

# --- GENERATE TABLE 2: Cluster vs. Healthy Controls ---
# Note: You must run the gsea_hc_list loop or topTable on fit_hc2 first to populate this
all_markers_hc <- lapply(1:4, function(i) {
  topTable(fit_hc2, coef = i, number = Inf)
})
names(all_markers_hc) <- c("1", "2", "3", "4")

table_hc_comparison <- process_marker_list(all_markers_hc)
write.csv(table_hc_comparison, "Table2_Cluster_vs_HealthyControls_Full.csv", row.names = FALSE)
#######################