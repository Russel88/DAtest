#' Log t-test

#' @export

DA.ltt <- function(count_table, outcome, delta = 1, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, paired = NULL, p.adj){
  
  tt <- function(x){
    tryCatch(t.test(x ~ outcome)$p.value, error = function(e){NA}) 
  }
  
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    outcome <- outcome[order(paired)]
    testStat <- testStat.pair
    tt <- function(x){
      tryCatch(t.test(x ~ outcome, paired = TRUE)$p.value, error = function(e){NA}) 
    }
  }
  
  count_table <- log(count_table+delta)
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  
  res <- data.frame(pval = apply(count.rel,1,tt))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  # Teststat
  outcome.num <- as.numeric(as.factor(outcome))-1
  testfun <- function(x){
    case <- x[outcome.num==1]
    control <- x[outcome.num==0]
    testStat(case,control) 
  }
  res$FC <- apply(count.rel,1,testfun)
  
  res$Feature <- rownames(res)
  res$Method <- "Log t-test"
  return(res)
}


