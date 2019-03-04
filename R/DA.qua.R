#' Quade test
#' 
#' Apply \code{quade.test} to multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param allResults If TRUE will return raw results from the \code{quade.test} function
#' @param ... Additional arguments for the \code{quade.test} function
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table, predictor, and paired variable
#' set.seed(4)
#' mat <- matrix(rnbinom(1500, size = 0.1, mu = 500), nrow = 100, ncol = 15)
#' rownames(mat) <- 1:100
#' pred <- c(rep("A", 5), rep("B", 5), rep("C", 5))
#' subject <- rep(1:5, 3)
#' 
#' # Running Quade test on each feature
#' res <- DA.qua(data = mat, predictor = pred, paired = subject)
#' @export

DA.qua <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  if(is.null(paired)) stop("Quade test needs a paired argument")
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
  } else {
    count_table <- data
  }
  
  # Define function
  qua <- function(x){
    tryCatch(quade.test(as.numeric(x), predictor, paired, ...), error = function(e){NA}) 
  }

  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  # Run tests
  reslist <- apply(count.rel,1,qua)

  # Extract results
  if(allResults){
    return(reslist)
  } else {
    res <- data.frame(statistic = sapply(reslist, function(x) x$statistic),
                      num.df = sapply(reslist, function(x) x$parameter[1]),
                      denom.df = sapply(reslist, function(x) x$parameter[2]),
                      pval = sapply(reslist, function(x) x$p.value))
    res[is.nan(res$statistic),"pval"] <- 1
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
    res$Feature <- gsub(".Quade.*","",rownames(res))
    res$Method <- "Quade (qua)" 
    
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
  
}
