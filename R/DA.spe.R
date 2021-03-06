#' Spearman's Rank Correlation
#' 
#' Apply Spearman correlation to multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param ... Additional arguments for the \code{cor.test} function
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 50, ncol = 20)
#' rownames(mat) <- 1:50
#' pred <- rnorm(20)
#' 
#' # Running Spearman correlation on each feature
#' res <- DA.spe(data = mat, predictor = pred)
#' @export

DA.spe <- function(data, predictor, relative = TRUE, p.adj = "fdr", ...){
  
  # Extract from phyloseq
  if(is(data, "phyloseq")){
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
  spe <- function(x){
    tryCatch(cor.test(x,predictor, method = "spearman", ...), error = function(e){NA}) 
  }
  
  # Run tests
  spes <- apply(count.rel, 1, spe)
  
  # Collect results
  res <- data.frame(pval = sapply(spes, function(x) x$p.value))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$rho <- sapply(spes, function(x) x$estimate)
  res$Feature <- rownames(res)
  res$Method <- "Spearman (spe)"
  
  if(is(data, "phyloseq")) res <- addTax(data, res)

  return(res)  
}

