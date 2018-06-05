#' Add tax_table from phyloseq object
#'
#' Internal function
#' @export

add.tax.DA <- function(data, res){
  
  if(!is.null(tax_table(data, errorIfNULL = FALSE))){
    tax <- unclass(tax_table(data))
    res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
    rownames(res) <- NULL
  } 
  
  return(res)

}

