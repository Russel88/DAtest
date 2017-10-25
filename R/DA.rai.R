#' RAIDA
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the raida function
#' @export

DA.rai <- function(data, predictor, p.adj = "fdr", ...){
  
  suppressMessages(library(RAIDA))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- unlist(sample_data(data)[,predictor])
  } else {
    count_table <- data
  }
  
  count_table.o <- as.data.frame(count_table[,order(predictor)])
  
  res <- raida(count_table.o, n.lib = as.numeric(table(predictor)), mtcm = p.adj, ...)
  res$Feature <- rownames(res)
  colnames(res)[1] <- "pval"
  colnames(res)[2] <- "pval.adj"
  res$Method <- "RAIDA (rai)" 
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  
  return(res)

}




