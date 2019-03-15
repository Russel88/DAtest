#' Friedman Rank Sum test
#' 
#' Apply friedman test to multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param allResults If TRUE will return raw results from the \code{friedman.test} function
#' @param ... Additional arguments for the \code{friedman.test} function
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table, predictor, and paired variable
#' set.seed(4)
#' mat <- matrix(rnbinom(1500, size = 0.1, mu = 500), nrow = 100, ncol = 15)
#' rownames(mat) <- 1:100
#' pred <- c(rep("A", 5), rep("B", 5), rep("C", 5))
#' subject <- rep(1:5, 3)
#' 
#' # Running Friedman test on each feature
#' res <- DA.fri(data = mat, predictor = pred, paired = subject)
#' @export

DA.fri <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  if(is.null(paired)) stop("Friedman test needs a paired argument")
  
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
  fri <- function(x){
    tryCatch(friedman.test(as.numeric(x), predictor, paired, ...), error = function(e){NA}) 
  }

  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  # Run tests
  reslist <- apply(count.rel,1,fri)

  # Collect results
  if(allResults){
    return(reslist)
  } else {
    res <- data.frame(statistic = sapply(reslist, function(x) x$statistic),
                      parameter = sapply(reslist, function(x) x$parameter),
                      pval = sapply(reslist, function(x) x$p.value))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
    res$Feature <- gsub(".Friedman.*","",rownames(res))
    res$Method <- "Friedman (fri)" 
    
    if(class(data) == "phyloseq") res <- addTax(data, res)
    
    return(res)
  }
  
}
