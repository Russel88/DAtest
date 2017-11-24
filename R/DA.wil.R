#' Wilcoxon Rank Sum and Signed Rank Test
#' 
#' Apply wilcoxon test for multiple features with one \code{predictor}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param testStat Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: \code{log2((mean(case abundances)+0.001)/(mean(control abundances)+0.001))}
#' @param testStat.pair Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: \code{log2(mean((case abundances+0.001)/(control abundances+0.001)))}
#' @param allResults If TRUE will return raw results from the \code{wilcox.test} function
#' @param ... Additional arguments for the \code{wilcox.test} function
#' @export

DA.wil <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", testStat = function(case,control){log2((mean(case)+0.001)/(mean(control)+0.001))}, testStat.pair = function(case,control){log2(mean((case+0.001)/(control+0.001)))}, allResults = FALSE, ...){
 
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
  wil <- function(x){
    tryCatch(wilcox.test(x ~ predictor, ...)$p.value, error = function(e){NA}) 
  }

  # Order data and define function for paired analysis
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    predictor <- predictor[order(paired)]
    testStat <- testStat.pair
    wil <- function(x){
      tryCatch(wilcox.test(x ~ predictor, paired = TRUE, ...)$p.value, error = function(e){NA}) 
    }
  }
  
  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  # Run tests
  if(allResults){
    if(is.null(paired)){
      wil <- function(x){
        tryCatch(wilcox.test(x ~ predictor, ...), error = function(e){NA}) 
      }
    } else {
      wil <- function(x){
        tryCatch(wilcox.test(x ~ predictor, paired = TRUE, ...), error = function(e){NA}) 
      }
    }
    return(apply(count.rel,1,wil))
  } else {
    res <- data.frame(pval = apply(count.rel,1,wil))
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
    res$Method <- "Wilcox (wil)" 
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }
 
}
