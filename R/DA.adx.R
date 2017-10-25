#' Aldex t.test and wilcox
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the aldex function
#' @export

DA.adx <- function(data, predictor, p.adj = "fdr", ...){
  
  suppressMessages(library(ALDEx2))
  
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
  
  x <- aldex(data.frame(count_table), predictor, verbose = FALSE, ...)
  x$we.ep.adj <- p.adjust(x$we.ep, method = p.adj)
  x$wi.ep.adj <- p.adjust(x$wi.ep, method = p.adj)
  x$Feature <- rownames(x)
  
  if(class(data) == "phyloseq") x <- add.tax.DA(data, x)
  
  return(x)
  
}




