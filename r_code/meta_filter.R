filter_by_metadata <- function(gse_id, expr_mat, pval_mat, threshold = 0.5) {
  message(paste("--- Processing", gse_id, "---"))
  
  # Load metadata
  gse_meta <- getGEO(gse_id, destdir = "transcriptome_data/uncompressed", getGPL = FALSE)
  metadata <- pData(gse_meta[[1]])
  
  # --- THE KEY FIX ---
  # Change row names from GSM IDs to the titles (3004A, etc.)
  rownames(metadata) <- metadata$title 
  
  # Flexible regex to catch "category: Kawasaki Disease"
  keep_idx <- which(grepl("Kawasaki|KD|Control|Healthy", 
                          metadata$characteristics_ch1, ignore.case = TRUE))
  
  # Subset matrices using the new row names
  meta_sub <- metadata[keep_idx, ]
  common_samples <- intersect(colnames(expr_mat), rownames(meta_sub))
  expr_sub <- expr_mat[, common_samples]
  pval_sub <- pval_mat[, common_samples]
  
  # 4. P-value Filtering
  # rowMeans now operates on a non-empty matrix
  keep_probes <- rowMeans(pval_sub < 0.05, na.rm = TRUE) > threshold
  
  if (sum(keep_probes) == 0) {
    warning(paste("No probes passed the p-value threshold for", gse_id))
    return(list(expr = expr_sub, metadata = meta_sub))
  }
  
  message(paste(gse_id, ": Found", length(common_samples), "samples. Kept", sum(keep_probes), "probes."))
  
  return(list(
    expr = expr_sub[keep_probes, ],
    metadata = meta_sub
  ))
}