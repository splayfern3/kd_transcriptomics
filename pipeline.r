library(tidyverse)
library(limma)
library(sva)
library(umap)
library(cluster)
library(msigdbr)
library(fgsea)

GSE73461_raw <- read.table(
  "transcriptome_data/uncompressed/GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cn <- colnames(GSE73461_raw)
is_id_col <- (cn == "ARRAY_ID")
is_pval_col <- grepl("Detection", cn)
is_expr_col <- !is_id_col & !is_pval_col

probe_ids <- GSE73461_raw$ARRAY_ID
expr_df <- GSE73461_raw[, is_expr_col, drop = FALSE]
pval_df <- GSE73461_raw[, is_pval_col, drop = FALSE]
GSE73461_expr <- as.matrix(expr_df)
GSE73461_pval <- as.matrix(pval_df)
rownames(GSE73461_expr) <- probe_ids
rownames(GSE73461_pval) <- probe_ids

pval_names_stripped <- sub("_Detection_Pval$", "", colnames(GSE73461_pval))
stopifnot(identical(colnames(GSE73461_expr), pval_names_stripped))
dim(GSE73461_expr)
dim(GSE73461_pval)
range(GSE73461_expr, na.rm = TRUE)
range(GSE73461_pval, na.rm = TRUE)

