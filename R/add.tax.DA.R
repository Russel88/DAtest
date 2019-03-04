#' Add tax_table from phyloseq object
#'
#' Internal function
#' @param data phyloseq object
#' @param res data.frame with results
#' @return \code{res} merged with the tax.table from the phyloseq object
#' @export
add.tax.DA <- function(data, res){
  
  loadNamespace("phyloseq")
  
  if(!is.null(phyloseq::tax_table(data, errorIfNULL = FALSE))){
    tax <- unclass(phyloseq::tax_table(data))
    res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
    rownames(res) <- NULL
  } 
  
  return(res)

}

