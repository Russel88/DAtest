#' EdgeR exact test
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the exactTest function
#' @export

DA.ere <- function(data, predictor, p.adj = "fdr", ...){
  
 library(edgeR)
  
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
  
  otu_table <- as.data.frame(count_table)
  x <- DGEList(counts = count_table, group = predictor, genes = data.frame(Feature = row.names(count_table)))
  x <- edgeR::calcNormFactors(x)
  x <- estimateTagwiseDisp(estimateCommonDisp(x))
  ta <- exactTest(x, ...)[[1]]
  colnames(ta)[3] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR exact"
 
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      ta <- merge(ta, tax, by.x = "Feature", by.y = "row.names")
      rownames(ta) <- NULL
    } 
  }
  
  return(ta) 
}


