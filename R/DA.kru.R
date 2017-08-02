#' Kruskal-Wallis test

#' @export

DA.kru <- function(count_table, outcome, p.adj, relative = TRUE){
 
  kru <- function(x){
    tryCatch(kruskal.test(as.numeric(x) ~ as.factor(outcome))$p.value, error = function(e){NA}) 
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
