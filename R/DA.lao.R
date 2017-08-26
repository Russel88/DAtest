#' ANOVA
#' 
#' With log transformation of counts before normalization.
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 1
#' @param ... Additional arguments for the aov functions
#' @export

DA.lao <- function(count_table, outcome,relative = TRUE, p.adj = "fdr", delta = 1, ...){
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(x ~ outcome, ...))[[1]][1,5]), error = function(e){NA}) 
  }
  
  count_table <- log(count_table+delta)
  if(relative) count_table <- apply(count_table,2,function(x) x/sum(x))
  
  res <- data.frame(pval = apply(count_table,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Log ANOVA"
  return(res)
}


