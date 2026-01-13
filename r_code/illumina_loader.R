load_illumina_geotxt <- function(path, floor_negatives_to_zero = FALSE) {
  # 1. Read the full file
  full_data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 2. SET ROW NAMES FIRST
  # Use the first column (ARRAY_ID) as row names and then remove it from the data frame
  # to prevent it from being accidentally treated as a numeric sample column.
  rownames(full_data) <- as.character(full_data[[1]])
  full_data <- full_data[, -1] 
  
  # 3. Identify P-value columns
  pval_cols <- grep("Pval|Detection", colnames(full_data), ignore.case = TRUE, value = TRUE)
  
  # 4. Identify Signal columns
  # Since we removed column 1, all remaining non-Pval columns are potential signals
  signal_cols <- gsub("_Detection_Pval|_Detection Pval|Detection Pval", "", pval_cols, ignore.case = TRUE)
  
  # 5. Extract and Convert to Matrix
  # Because we set rownames in step 2, expr and pval will inherit them automatically
  expr <- as.matrix(full_data[, signal_cols])
  pval <- as.matrix(full_data[, pval_cols])
  
  # 6. Optional negative flooring
  if(floor_negatives_to_zero) {
    expr[expr < 0] <- 0
  }
  
  # Consistency Check
  if(ncol(expr) != ncol(pval)) {
    stop(paste("Mismatch! Found", ncol(expr), "signals and", ncol(pval), "p-values."))
  }
  
  message("Successfully loaded ", ncol(expr), " samples and ", nrow(expr), " probes.")
  
  return(list(expr_lin = expr, pval = pval))
}