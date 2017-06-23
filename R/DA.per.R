#' Permutation test
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' Pvalues are now two-sided, and test statistic is log fold change

#' @export

DA.per <- function(otu_table, outcome, noOfIterations = 10000, seed = as.numeric(Sys.time()), margin = 50, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}){

  outcome <- as.numeric(as.factor(outcome))-1
  
  otu_table <- apply(otu_table, 2, function(x) x/sum(x))
  
  set.seed(seed)
  nullStatList <- list()
  
  # Create shuffled outcomes
  shuffledOutcomesList <- list()
  for (k in 1:noOfIterations){
    shuffledOutcomesList[[k]] <- sample(outcome)
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
        ptemp <- mean(abs(nullStatTemp) >= abs(realStat))
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
    
    p[i]         <- mean(abs(nullStat) >= abs(realStat))
    if(p[i] == 0){ #only at level of support
      p[i] <- 1/noOfIterations
    }
    FC[i]        <- realStat
  }
  output_df <- data.frame(OTU = row.names(otu_table), pval = p, FC, coverage)
  output_df$Method <- "Permutation"
  return(output_df)
}


