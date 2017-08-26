#' Welch t-test
#
# With log transformaiton of counts before normalization
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for log transformation.
#' @param testStat Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param testStat.pair Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: mean(log((case abundances+1)/(control abundances+1)))
#' @param ... Additional arguments for the t.test function
#' @export

DA.ltt <- function(count_table, outcome, paired = NULL, relative = TRUE, p.adj = "fdr", delta = 1, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, ...){
  
  tt <- function(x){
    tryCatch(t.test(x ~ outcome, ...)$p.value, error = function(e){NA}) 
  }
  
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    outcome <- outcome[order(paired)]
    testStat <- testStat.pair
    tt <- function(x){
      tryCatch(t.test(x ~ outcome, paired = TRUE, ...)$p.value, error = function(e){NA}) 
    }
  }
  
  count_table <- log(count_table+delta)
  if(relative) count_table <- apply(count_table,2,function(x) x/sum(x))
  
  res <- data.frame(pval = apply(count_table,1,tt))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  # Teststat
  outcome.num <- as.numeric(as.factor(outcome))-1
  testfun <- function(x){
    case <- x[outcome.num==1]
    control <- x[outcome.num==0]
    testStat(case,control) 
  }
  res$FC <- apply(count_table,1,testfun)
  
  res$Feature <- rownames(res)
  res$Method <- "Log t-test"
  return(res)
}


