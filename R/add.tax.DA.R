addTax <- function(data, res){
  
  loadNamespace("phyloseq")
  
  if(!is.null(phyloseq::tax_table(data, errorIfNULL = FALSE))){
    tax <- unclass(phyloseq::tax_table(data))
    res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
    rownames(res) <- NULL
  } 
  
  return(res)

}

