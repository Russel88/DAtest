#' EdgeR exact test - RLE normalization
#'
#' Implementation of edgeR exact test for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param ... Additional arguments for the \code{calcNormFactors}, \code{estimateCommonDisp}, \code{estimateTagwiseDisp} and \code{exactTest} functions
#' @return A data.frame with with results.
#' @export

DA.ere2 <- function(data, predictor, p.adj = "fdr", ...){
  
  ok <- tryCatch({
    loadNamespace("edgeR")
    TRUE
  }, error=function(...) FALSE)
  
  if (ok){
    # Extract from phyloseq
    if(class(data) == "phyloseq"){
      DAdata <- DA.phyloseq(data, predictor)
      count_table <- DAdata$count_table
      predictor <- DAdata$predictor
    } else {
      count_table <- data
    }
    
    count_table <- as.data.frame(count_table)
    x <- edgeR::DGEList(counts = count_table, group = predictor, genes = data.frame(Feature = row.names(count_table)))
    
    # Extract arguments
    DA.ere.args <- list(...)
    calcNormFactors.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(edgeR::calcNormFactors))]
    estimateCommonDisp.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(edgeR::estimateCommonDisp))]
    estimateTagwiseDisp.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(edgeR::estimateTagwiseDisp))]
    exactTest.args <- DA.ere.args[names(DA.ere.args) %in% names(formals(edgeR::exactTest))]
    
    # Normalize and fit model
    x <- do.call(edgeR::calcNormFactors,c(list(x, method = "RLE"),calcNormFactors.args))
    x <- do.call(edgeR::estimateCommonDisp,c(list(x),estimateCommonDisp.args))
    x <- do.call(edgeR::estimateTagwiseDisp,c(list(x),estimateTagwiseDisp.args))
    ta <- do.call(edgeR::exactTest,c(list(x),exactTest.args))[[1]]
    
    colnames(ta)[3] <- "pval"
    ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
    ta$ordering <- NA
    ta[!is.na(ta$logFC) & ta$logFC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
    ta[!is.na(ta$logFC) & ta$logFC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
    ta$Feature <- rownames(ta)
    ta$Method <- "EdgeR exact - RLE (ere2)"
    
    if(class(data) == "phyloseq") ta <- add.tax.DA(data, ta)
    
    return(ta) 
  } else {
    stop("edgeR package required")
  }
  

}


