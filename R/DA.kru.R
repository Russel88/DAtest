#' Kruskal-Wallis test
#' 
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the kruskal.test function
#' @export

DA.kru <- function(count_table, outcome, relative = TRUE, p.adj = "fdr", ...){
 
  kru <- function(x){
    tryCatch(kruskal.test(as.numeric(x) ~ as.factor(outcome), ...)$p.value, error = function(e){NA}) 
  }

  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  res <- data.frame(pval = apply(count.rel,1,kru))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Kruskal-Wallis" 
  return(res)
 
}
