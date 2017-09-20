#' Spearman's Rank Correlation
#'
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the cor.test function
#' @export

DA.spe <- function(data, predictor, relative = TRUE, p.adj = "fdr", ...){
  
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
  
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  spe <- function(x){
    tryCatch(cor.test(x,predictor, method = "spearman", ...)$p.value, error = function(e){NA}) 
  }
  
  spe.cor <- function(x){
    tryCatch(cor(x,predictor, method = "spearman"), error = function(e){NA}) 
  }

  res <- data.frame(pval = apply(count.rel,1,spe))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$rho <- apply(count.rel,1,spe.cor)
  res$Feature <- rownames(res)
  res$Method <- "Spearman (spe)"
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)

  return(res)  
}

