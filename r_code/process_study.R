# --- DATA PRE-PROCESSING: FILTERING & COLLAPSING ---

# 1. Define a robust helper function for this specific study
# This ensures we handle the "category: KD" and "category: Control" labeling correctly
process_study <- function(expr_mat, metadata, mapped_df, dataset_name) {
  message(paste("Processing", dataset_name, "..."))
  
  # Identify target columns (KD and Control only)
  # We look for 'KD' or 'Kawasaki' and 'Control' or 'Healthy'
  keep_idx <- which(
    grepl("KD\\b|Kawasaki", metadata$characteristics_ch1, ignore.case = TRUE) | 
      grepl("Control|Healthy", metadata$characteristics_ch1, ignore.case = TRUE)
  )
  
  meta_sub <- metadata[keep_idx, ]
  
  # Subset matrix columns to match clean metadata
  expr_sub <- expr_mat[, colnames(expr_mat) %in% meta_sub$title]
  
  # Align Probes (Annotation)
  mapped_sub <- mapped_df[mapped_df$PROBEID %in% rownames(expr_sub), ]
  expr_ready <- expr_sub[as.character(mapped_sub$PROBEID), ]
  
  # Collapse Rows to Gene Symbols using MaxMean
  collapsed <- WGCNA::collapseRows(
    datET = expr_ready,
    rowGroup = mapped_sub$ENTREZID,
    rowID = mapped_sub$PROBEID,
    method = "MaxMean"
  )$datETcollapsed
  
  message(paste(dataset_name, "Complete. Samples:", ncol(collapsed), "Genes:", nrow(collapsed)))
  return(collapsed)
}
