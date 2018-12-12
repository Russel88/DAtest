#' Welch t-test - Multiplicative zero-correction and center log-ratio normalization. 
#' 
#' Apply welch t-test to multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param delta Numeric. Pseudocount for zero-correction.
#' @param testStat Function. Function for calculating fold change. Should take two vectors as arguments. Default is a simple difference: \code{mean(case abundances)-mean(control abundances)}
#' @param testStat.pair Function. Function for calculating fold change. Should take two vectors as arguments. Default is a simple difference: \code{mean(case abundances-control abundances)}
#' @param allResults If TRUE will return raw results from the \code{t.test} function
#' @param ... Additional arguments for the \code{t.test} function
#' @export

DA.ttc <- function(data, predictor, paired = NULL, p.adj = "fdr", delta = 1, testStat = function(case,control){mean(case)-mean(control)}, testStat.pair = function(case,control){mean(case-control)}, allResults = FALSE, ...){
  
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
  tt <- function(x){
    tryCatch(t.test(x ~ predictor, ...)$p.value, error = function(e){NA}) 
  }

  # Order data and define function for paired analysis
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    predictor <- predictor[order(paired)]
    testStat <- testStat.pair
    tt <- function(x){
      tryCatch(t.test(x ~ predictor, paired = TRUE, ...)$p.value, error = function(e){NA}) 
    }
  }
  
  # Zero-correction
  count_table <- apply(count_table, 2, function(y) sapply(y,function(x) ifelse(x==0,delta,(1-(sum(y==0)*delta)/sum(y))*x)))
  if(any(count_table <= 0)) stop("count_table should only contain positive values")
  
  # CLR transformation
  count_table <- norm_clr(count_table)
  
  # Run tests
  if(allResults){
    if(is.null(paired)){
      tt <- function(x){
        tryCatch(t.test(x ~ predictor, ...), error = function(e){NA}) 
      }
    } else {
      tt <- function(x){
        tryCatch(t.test(x ~ predictor, paired = TRUE, ...), error = function(e){NA}) 
      }
    }
    return(apply(count_table,1,tt))
  } else {
    res <- data.frame(pval = apply(count_table,1,tt))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    # Teststat
    predictor.num <- as.numeric(as.factor(predictor))-1
    testfun <- function(x){
      case <- x[predictor.num==1]
      control <- x[predictor.num==0]
      testStat(case,control) 
    }
    res$log2FC <- apply(count_table,1,testfun)
    res$ordering <- NA
    res[!is.na(res$log2FC) & res$log2FC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
    res[!is.na(res$log2FC) & res$log2FC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
    res$Feature <- rownames(res)
    res$Method <- "t-test - CLR (ttc)"
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
}
