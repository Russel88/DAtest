#' Log ANOVA

#' @export

DA.lao <- function(count_table, outcome, delta = 1, p.adj){
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(x ~ outcome))[[1]][1,5]), error = function(e){NA}) 
  }
  
  count_table <- log(count_table+delta)
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  
  res <- data.frame(pval = apply(count.rel,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Log ANOVA"
  return(res)
}


