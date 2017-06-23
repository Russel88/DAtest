#' Log t-test2

#' @export

DA.ltt2 <- function(otu_table, outcome, delta = 0.001, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}){
  
  tt <- function(x){
    tryCatch(t.test(x ~ outcome)$p.value, error = function(e){NA}) 
  }
  otu.rel <- apply(otu_table,2,function(x) x/sum(x))
  otu.rel <- log(otu.rel+delta)
  res <- data.frame(pval = apply(otu.rel,1,tt))
  
  # Teststat
  outcome.num <- as.numeric(as.factor(outcome))-1
  testfun <- function(x){
    case <- x[outcome.num==1]
    control <- x[outcome.num==0]
    testStat(case,control) 
  }
  res$FC <- apply(otu.rel,1,testfun)
  
  res$OTU <- rownames(res)
  res$Method <- "Log t-test2"
  return(res)
}
