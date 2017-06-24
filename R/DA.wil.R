#' Wilcox test

#' @export

DA.wil <- function(otu_table, outcome, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, p.adj){
 
  wil <- function(x){
    tryCatch(wilcox.test(x ~ outcome)$p.value, error = function(e){NA}) 
  }

  otu.rel <- apply(otu_table,2,function(x) x/sum(x))
  res <- data.frame(pval = apply(otu.rel,1,wil))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  # Teststat
  outcome.num <- as.numeric(as.factor(outcome))-1
  testfun <- function(x){
    case <- x[outcome.num==1]
    control <- x[outcome.num==0]
    testStat(case,control) 
  }
  
  res$FC <- apply(otu.rel,1,testfun)
  
  res$OTU <- rownames(res)
  res$Method <- "Wilcox" 
  return(res)
 
}
