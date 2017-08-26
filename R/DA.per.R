#' Permutation test of user-defined test statistic
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' P-values are now two-sided, and test statistic is a simple log fold change
#' 
#' A paired permutation test is implemented specifically for this package. The test is similar to the original, but with a different test statistic and permutation scheme. The permutations are constrained in the paired version such that the outcome is only permuted within each level of the paired argument (e.g. subjects). The test statistic first finds the log-ratio between the two outcome levels (e.g. case and control) for each level of the paired argument and the final statistic is the mean of these log-ratios.
#' 
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param testStat Function. Function for the test statistic. Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param testStat.pair Function. Function for test statistc for paired analysis. Should take two vectors as arguments. Default is a log fold change: mean(log((case abundances+1)/(control abundances+1)))
#' @param noOfIterations Integer. Iterations for permutations. Default 10000
#' @param margin Numeric. Margin for when to stop iterations if p-value is high and unlikely to become low
#' @export

DA.per <- function(count_table, outcome, paired = NULL, relative = TRUE, p.adj = "fdr", testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, noOfIterations = 10000, margin = 50){

  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    outcome <- outcome[order(paired)]
    testStat <- testStat.pair
  }
  
  outcome <- as.numeric(as.factor(outcome))-1
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
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
  
  iterations <- nrow(count_table)
  p <- numeric(iterations)
  FC <- numeric(iterations)
  coverage <- numeric(iterations)
  
  for(i in 1:iterations){
    
    count_row      <- as.numeric(count_table[i,])
    
    real_case    <- count_row[outcome==1]
    real_control <- count_row[outcome==0]
    realStat     <- testStat(real_case,real_control) 
    
    Wnull <- numeric(noOfIterations)
    
    for(j in 1:noOfIterations){
      case     <- count_row[shuffledOutcomesList[[j]]==1]
      control  <- count_row[shuffledOutcomesList[[j]]==0]
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
  output_df <- data.frame(Feature = row.names(count_table), pval = p, FC, coverage)
  output_df$pval.adj <- p.adjust(output_df$pval, method = p.adj)
  output_df$Method <- "Permutation"
  return(output_df)
}


