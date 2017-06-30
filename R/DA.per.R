#' Permutation test
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' 
#' Pvalues are now two-sided, and test statistic is a simple log fold change
#' 
#' A paired permutaition test is implemented specifically for this package. The test is similar to the original, but with a different test statistic and permutaition scheme. The permutations are constrained in the paired version such that the outcome is only permuted within each level of the paired argument (e.g. subjects). The test statistic first finds the log-ratio between the two outcome levels (e.g. case and control) for each level of the paired argument and the final statistic is the mean of these log-ratios.

#' @export

DA.per <- function(otu_table, outcome, paired = NULL, noOfIterations = 10000, seed = as.numeric(Sys.time()), margin = 50, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, p.adj){

  if(!is.null(paired)){
    otu_table <- otu_table[,order(paired)]
    outcome <- outcome[order(paired)]
    testStat <- testStat.pair
  }
  
  outcome <- as.numeric(as.factor(outcome))-1
  otu_table <- apply(otu_table, 2, function(x) x/sum(x))
  
  set.seed(seed)
  nullStatList <- list()
  
  # Create shuffled outcomes
  shuffledOutcomesList <- list()
  if(is.null(paired)){
    for (k in 1:noOfIterations){
      shuffledOutcomesList[[k]] <- sample(outcome)
    }
  } else {
    for (k in 1:noOfIterations){
      shuffledOutcomesList[[k]] <- unlist(lapply(1:(length(outcome)/2),function(x) sample(c(0,1))))
    }
  }
  
  iterations <- nrow(otu_table)
  p <- numeric(iterations)
  FC <- numeric(iterations)
  coverage <- numeric(iterations)
  
  for(i in 1:iterations){
    
    otu_row      <- as.numeric(otu_table[i,])
    
    real_case    <- otu_row[outcome==1]
    real_control <- otu_row[outcome==0]
    realStat     <- testStat(real_case,real_control) 
    
    Wnull <- numeric(noOfIterations)
    
    for(j in 1:noOfIterations){
      case     <- otu_row[shuffledOutcomesList[[j]]==1]
      control  <- otu_row[shuffledOutcomesList[[j]]==0]
      Wnull[j] <- testStat(case,control)
      
      if(j %in% 10^(1:100)){
        nullStatTemp <- Wnull[1:j]
        ptemp <- sum(abs(nullStatTemp) >= abs(realStat))/(noOfIterations+1)
        if(ptemp > (margin*(1/j))){
          nullStat <- Wnull[1:j]
          coverage[i] <- j
          break
        }
      }
      if(j == noOfIterations){
        coverage[i] <- j
        nullStat <- Wnull
      }
    }
    
    p[i]         <- sum(abs(nullStat) >= abs(realStat))/(noOfIterations+1)
    if(p[i] == 0){ #only at level of support
      p[i] <- 1/noOfIterations
    }
    FC[i]        <- realStat
  }
  output_df <- data.frame(OTU = row.names(otu_table), pval = p, FC, coverage)
  output_df$pval.adj <- p.adjust(output_df$pval, method = p.adj)
  output_df$Method <- "Permutation"
  return(output_df)
}


