#' Permutation test of user-defined test statistic
#'
#' Test significance of multiple features with a binary \code{predictor} with a permutation scheme.
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' P-values are now two-sided, and test statistic is a simple log fold change
#' 
#' A paired permutation test is implemented specifically for this package. The test is similar to the original, but with a different test statistic and permutation scheme. The permutations are constrained in the paired version such that the predictor is only permuted within each level of the paired argument (e.g. subjects). The test statistic first finds the log-ratio between the two predictor levels (e.g. case and control) for each level of the paired argument and the final statistic is the mean of these log-ratios.
#' 
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param relative Logical. Should \code{data} be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param testStat Function. Function for the test statistic. Should take two vectors as arguments. Default is a log fold change: \code{log2((mean(case abundances)+1)/(mean(control abundances)+1))}
#' @param testStat.pair Function. Function for test statistic for paired analysis. Should take two vectors as arguments. Default is a log fold change: \code{log2(mean((case abundances+1)/(control abundances+1)))}
#' @param noOfIterations Integer. Iterations for permutations. Default 10000
#' @param margin Numeric. Margin for when to stop iterations if p-value is high and unlikely to become low
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 100, ncol = 10)
#' rownames(mat) <- 1:100
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running permutation test on each feature
#' res <- DA.per(data = mat, predictor = pred)
#' @export

DA.per <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", testStat = function(case,control){log2((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){log2(mean((case+1)/(control+1)))}, noOfIterations = 10000, margin = 50){

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
  } else {
    count_table <- data
  }

  # Order data for paired analysis
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    predictor <- predictor[order(paired)]
    testStat <- testStat.pair
  }
  
  pred.lev <- levels(as.factor(predictor)) 
  predictor <- as.numeric(as.factor(predictor))-1
  
  # Relative abundance
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  nullStatList <- list()
  
  # Create shuffled predictors
  shuffledpredictorsList <- list()
  if(is.null(paired)){
    for (k in seq_len(noOfIterations)){
      shuffledpredictorsList[[k]] <- sample(predictor)
    }
  } else {
    for (k in seq_len(noOfIterations)){
      shuffledpredictorsList[[k]] <- unlist(lapply(seq_len(length(predictor)/2),function(x) sample(c(0,1))))
    }
  }
  
  iterations <- nrow(count.rel)
  p <- numeric(iterations)
  FC <- numeric(iterations)
  coverage <- numeric(iterations)
  
  # For each feature
  for(i in seq_len(iterations)){
    
    count_row      <- as.numeric(count.rel[i,])
    
    real_case    <- count_row[predictor==1]
    real_control <- count_row[predictor==0]
    realStat     <- testStat(real_case,real_control) 
    
    Wnull <- numeric(noOfIterations)
    
    # For each iteration
    for(j in seq_len(noOfIterations)){
      case     <- count_row[shuffledpredictorsList[[j]]==1]
      control  <- count_row[shuffledpredictorsList[[j]]==0]
      Wnull[j] <- testStat(case,control)
      
      # Break if p-value is high (as defined by margin)
      if(j %in% 10^(1:100)){
        nullStatTemp <- Wnull[seq_len(j)]
        ptemp <- sum(abs(nullStatTemp) >= abs(realStat))/(noOfIterations+1)
        if(ptemp > (margin*(1/j))){
          nullStat <- Wnull[seq_len(j)]
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
  
  # Collect results
  output_df <- data.frame(Feature = row.names(count_table), pval = p, log2FC = FC, coverage)
  output_df$ordering <- NA
  output_df[!is.na(output_df$log2FC) & output_df$log2FC > 0,"ordering"] <- paste0(pred.lev[2],">",pred.lev[1])
  output_df[!is.na(output_df$log2FC) & output_df$log2FC < 0,"ordering"] <- paste0(pred.lev[1],">",pred.lev[2])
  output_df$pval.adj <- p.adjust(output_df$pval, method = p.adj)
  output_df$Method <- "Permutation (per)"
  
  if(class(data) == "phyloseq") output_df <- add.tax.DA(data, output_df)
  
  return(output_df)
}


