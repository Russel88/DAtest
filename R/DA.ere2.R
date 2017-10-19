#' EdgeR exact test - RLE normalization
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the calcNormFactors, estimateCommonDisp, estimateTagwiseDisp and exactTest functions
#' @export

DA.ere2 <- function(data, predictor, p.adj = "fdr", ...){
  
  suppressMessages(library(edgeR))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
  } else {
    count_table <- data
  }
  
  count_table <- as.data.frame(count_table)
  x <- DGEList(counts = count_table, group = predictor, genes = data.frame(Feature = row.names(count_table)))
  
  DA.ere.args <- list(...)
  calcNormFactors.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(calcNormFactors))]
  estimateCommonDisp.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(estimateCommonDisp))]
  estimateTagwiseDisp.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(estimateTagwiseDisp))]
  exactTest.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(exactTest))]
  
  x <- do.call(edgeR::calcNormFactors,c(list(x, method = "RLE"),calcNormFactors.args))
  x <- do.call(estimateCommonDisp,c(list(x),estimateCommonDisp.args))
  x <- do.call(estimateTagwiseDisp,c(list(x),estimateTagwiseDisp.args))
  ta <- do.call(exactTest,c(list(x),exactTest.args))[[1]]
  colnames(ta)[3] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR exact - RLE (ere2)"
 
  if(class(data) == "phyloseq") ta <- add.tax.DA(data, ta)
  
  return(ta) 
}


