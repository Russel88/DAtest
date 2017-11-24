#' Pearson Correlation
#' 
#' Apply pearson correlation between multiple features and one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param ... Additional arguments for the \code{cor.test} function
#' @export

DA.pea <- function(data, predictor, relative = TRUE, p.adj = "fdr", ...){
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
  } else {
    count_table <- data
  }

  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  # Define function
  pea <- function(x){
    tryCatch(cor.test(x,predictor, ...), error = function(e){NA}) 
  }
  
  # Run tests
  peas <- apply(count.rel, 1, pea)
  
  # Collect results
  res <- data.frame(pval = sapply(peas, function(x) x$p.value))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$cor <- sapply(peas, function(x) x$estimate)
  res$`cor_2.5%` <- sapply(peas, function(x) x$conf.int[1])
  res$`cor_97.5%` <- sapply(peas, function(x) x$conf.int[2])
  res$Feature <- rownames(res)
  res$Method <- "Pearson (pea)"

  # Add tax table if available
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)

  return(res)  
}

