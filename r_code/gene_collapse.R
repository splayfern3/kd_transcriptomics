collapse <- function(expr, anno) {
  expr_s <- expr[as.character(anno$PROBEID)] # align expr to anno
  
  collapsed <- collapseRows(
    datET = expr_s,
    rowGroup = anno$ENTREZID,
    rowID = anno$PROBEID
    method = "MaxMean"
  )
  
  return(collapsed)
}