#' Spearman's Rank Correlation
#'
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Numeric. The outcome of interest. E.g. case and control
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @export

DA.spe <- function(count_table, outcome, relative = TRUE, p.adj = "fdr", ...){
  
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  spe <- function(x){
    tryCatch(cor.test(x,outcome, method = "spearman", ...)$p.value, error = function(e){NA}) 
  }
  
  spe.cor <- function(x){
    tryCatch(cor(x,outcome, method = "spearman"), error = function(e){NA}) 
  }

  res <- data.frame(pval = apply(count.rel,1,spe))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$rho <- apply(count.rel,1,spe.cor)
  res$Feature <- rownames(res)
  res$Method <- "Spearman"

  return(res)  
}

