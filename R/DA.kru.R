#' Kruskal-Wallis test
#' 
#' Apply kruskal-wallis test on multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a phyloseq object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param allResults If TRUE will return raw results from the \code{kruskal.test} function
#' @param ... Additional arguments for the \code{kruskal.test} function
#' @export

DA.kru <- function(data, predictor, relative = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
  } else {
    count_table <- data
  }

  # Define function
  kru <- function(x){
    tryCatch(kruskal.test(as.numeric(x) ~ predictor, ...)$p.value, error = function(e){NA}) 
  }

  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  # Run tests
  if(allResults){
    ## Define function for allResults 
    if(is.null(paired)){
      kru <- function(x){
        tryCatch(kruskal.test(x ~ predictor, ...), error = function(e){NA}) 
      }
    } else {
      kru <- function(x){
        tryCatch(kruskal.test(x ~ predictor, paired = TRUE, ...), error = function(e){NA}) 
      }
    }
    return(apply(count.rel,1,kru))
  } else {
    res <- data.frame(pval = apply(count.rel,1,kru))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- rownames(res)
    res$Method <- "Kruskal-Wallis (kru)" 
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
 
}
