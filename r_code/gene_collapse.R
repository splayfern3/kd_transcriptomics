collapse <- function(expr, anno) {
  # 1. Filter annotation to only include probes that actually exist in the expr matrix
  # This prevents 'subscript out of bounds' or empty matches
  anno_sub <- anno[anno$PROBEID %in% rownames(expr), ]
  
  # 2. Align expr to anno - THE COMMA IS CRITICAL HERE
  # expr[rows, columns]. Leaving columns blank selects all samples.
  expr_s <- as.matrix(expr[as.character(anno_sub$PROBEID), ])
  
  # 3. Run WGCNA collapse logic
  collapsed_obj <- WGCNA::collapseRows(
    datET = expr_s,
    rowGroup = anno_sub$SYMBOL,
    rowID = anno_sub$PROBEID,
    method = "MaxMean"
  )
  
  # 4. Return the matrix specifically so it can be piped
  return(collapsed_obj$datETcollapsed)
}