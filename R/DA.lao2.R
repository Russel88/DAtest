#' ANOVA
#' 
#' With log transformation of relative abundances.
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 0.001
#' @param ... Additional arguments for the aov functions
#' @export

DA.lao2 <- function(count_table, outcome, p.adj = "fdr", delta = 0.001, ...){
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(x ~ outcome, ...))[[1]][1,5]), error = function(e){NA}) 
  }
  
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  res <- data.frame(pval = apply(count.rel,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Log ANOVA 2"
  return(res)
}


