#' ANOVA

#' @export

DA.aov <- function(count_table, outcome, p.adj, relative){
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(x ~ outcome))[[1]][1,5]), error = function(e){NA}) 
  }
  
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  res <- data.frame(pval = apply(count.rel,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "ANOVA"
  return(res)
}


