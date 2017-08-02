#' Log ANOVA 2

#' @export

DA.lao2 <- function(count_table, outcome, delta = 0.001, p.adj){
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(x ~ outcome))[[1]][1,5]), error = function(e){NA}) 
  }
  
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  res <- data.frame(pval = apply(count.rel,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Log ANOVA 2"
  return(res)
}


