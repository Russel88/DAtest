#' ANOVA
#' 
#' Apply ANOVA to multiple features with one \code{predictor}, with log transformation of relative abundances.
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 0.001
#' @param allResults If TRUE will return raw results from the \code{aov} function
#' @param ... Additional arguments for the \code{aov} functions
#' @return A data.frame with with results.
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1500, size = 0.1, mu = 500), nrow = 100, ncol = 15)
#' rownames(mat) <- 1:100
#' pred <- c(rep("A", 5), rep("B", 5), rep("C", 5))
#' 
#' # Running ANOVA on each feature
#' res <- DA.lao2(data = mat, predictor = pred)
#' @export

DA.lao2 <- function(data, predictor, covars = NULL, p.adj = "fdr", delta = 0.001, allResults = FALSE, ...){
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired = NULL, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    covars <- DAdata$covars
  } else {
    count_table <- data
  }
  if(!is.null(covars)){
    for(i in seq_along(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }
  
  # Define model
  if(is.null(covars)){
    form <- paste("x ~ predictor")
  } else {
    form <- paste("x ~ ",paste(names(covars), collapse="+"),"+ predictor",sep = "")
  }
  
  # Define function
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(as.formula(form), ...))[[1]][(length(covars)+1),5]), error = function(e){NA}) 
  }
  
  # Relative abundance and log
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  # Run tests
  if(allResults){
    ao <- function(x){
      tryCatch(aov(as.formula(form), ...), error = function(e){NA}) 
    }
    return(apply(count.rel,1,ao))
  } else {
    res <- data.frame(pval = apply(count.rel,1,ao))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- rownames(res)
    res$Method <- "Log ANOVA 2 (lao2)"
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
}


