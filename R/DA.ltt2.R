#' Welch t-test
#' 
#' Apply welch t-test to multiple features and one \code{predictor}, and with log transformed relative abundances
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param delta Numeric. Pseudocount for log transformation. Default 0.001
#' @param testStat Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: \code{log2(mean(exp(case abundances)) / mean(exp(control abundances)))}
#' @param testStat.pair Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: \code{log2(mean(exp(case abundances) / exp(control abundances)))}
#' @param allResults If TRUE will return raw results from the \code{t.test} function
#' @param ... Additional arguments for the \code{t.test} function
#' @return A data.frame with with results.
#' @export

DA.ltt2 <- function(data, predictor, paired = NULL, p.adj = "fdr", delta = 0.001, testStat = function(case,control){log2((mean(exp(case)))/(mean(exp(control))))}, testStat.pair = function(case,control){log2(mean((exp(case))/(exp(control))))},allResults = FALSE, ...){
  
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
  
  # Relative abundance and log
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel+delta)
  
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
    return(apply(count.rel,1,tt))
  } else {
    res <- data.frame(pval = apply(count.rel,1,tt))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    # Teststat
    predictor.num <- as.numeric(as.factor(predictor))-1
    testfun <- function(x){
      case <- x[predictor.num==1]
      control <- x[predictor.num==0]
      testStat(case,control) 
    }
    res$log2FC <- apply(count.rel,1,testfun)
    res$ordering <- NA
    res[!is.na(res$log2FC) & res$log2FC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
    res[!is.na(res$log2FC) & res$log2FC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
    res$Feature <- rownames(res)
    res$Method <- "Log t-test2 (ltt2)"
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
}
