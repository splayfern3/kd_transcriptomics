# R/loaders_illumina.R

#Function to load illumina data
load_illumina_geotxt <- function(path, 
                                 id_col = "raw.ARRAY_ID",
                                 detect_pattern = "Detection",
                                 pval_suffix_regex = "_Detection_Pval$",
                                 floor_negatives_to_zero = TRUE) {
  
  raw <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cn <- colnames(raw)
  is_id_col <- (cn == id_col)
  is_pval_col <- grepl(detect_pattern, cn)
  is_expr_col <- !is_id_col & !is_pval_col
  
  probe_ids <- raw[[id_col]]
  
  expr <- as.matrix(raw[, is_expr_col, drop = FALSE]) #expression data in matrix format
  pval <- as.matrix(raw[, is_pval_col, drop = FALSE]) #pvalue data in matrix format
  #set probe ids to the rows
  rownames(expr) <- probe_ids 
  rownames(pval) <- probe_ids
  
  stopifnot(ncol(expr) == ncol(pval))
  
  pval_names_stripped <- sub(pval_suffix_regex, "", colnames(pval)) #rename pval data by removing the detection_pval part
  stopifnot(identical(colnames(expr), pval_names_stripped))
  
  if (floor_negatives_to_zero) {
    expr[expr < 0] <- 0
  }
  #returns
  list(
    raw = raw, #return raw data
    expr_lin = expr, #return expr matrix
    pval = pval, #return p value matrix
    path = path #return path of files
  )
}
