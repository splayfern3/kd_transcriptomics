filter_by_metadata <- function(gse_id, expr, pval, threshold = 0.5, data_dir = "transcriptome_data/uncompressed") {
  library(GEOquery)
  library(dplyr)
  
  message(paste("--- Processing", gse_id, "---"))
  
  # 1. Load Metadata
  gse_meta_list <- getGEO(gse_id, destdir = data_dir, getGPL = FALSE)
  metadata <- pData(gse_meta_list[[1]])
  rm(gse_meta_list); gc()
  
  # 2. Identify Groups
  metadata_clean <- dplyr::select(metadata, title, contains("characteristics_ch1"))
  colnames(metadata_clean)[2] <- "Diagnosis"
  if(ncol(metadata_clean) >= 3) colnames(metadata_clean)[3] <- "Stage"
  
  # 3. Target Indices
  target_indices <- which(
    grepl("Kawasaki", metadata_clean$Diagnosis, ignore.case = TRUE) | 
      grepl("Control", metadata_clean$Diagnosis, ignore.case = TRUE)
  )
  
  # --- DIAGNOSTIC CHECK ---
  message(paste("Target samples found:", length(target_indices)))
  
  # 4. P-value Thresholding
  pval_subset <- pval[, target_indices]
  
  # Check for NAs
  if(all(is.na(pval_subset))) {
    warning(paste("CRITICAL: All P-values are NA for", gse_id, ". Check your loader function! Returning unfiltered data."))
    return(list(expr = expr[, target_indices], metadata = metadata_clean[target_indices, ]))
  }
  
  is_detected <- pval_subset < 0.05
  is_detected[is.na(is_detected)] <- FALSE # Treat NAs as "Not Detected"
  
  sub_meta <- metadata_clean[target_indices, ]
  kd_idx <- which(grepl("Kawasaki", sub_meta$Diagnosis, ignore.case = TRUE))
  con_idx <- which(grepl("Control", sub_meta$Diagnosis, ignore.case = TRUE))
  
  # 50% Rule
  keep_kd <- rowSums(is_detected[, kd_idx, drop=FALSE]) >= (length(kd_idx) * threshold)
  keep_con <- rowSums(is_detected[, con_idx, drop=FALSE]) >= (length(con_idx) * threshold)
  keep_probes <- keep_kd | keep_con
  
  message(paste(gse_id, ": Kept", sum(keep_probes), "probes out of", nrow(expr)))
  
  return(list(expr = expr[keep_probes, target_indices], metadata = sub_meta))
}