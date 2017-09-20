#' Permutation test of user-defined test statistic
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' P-values are now two-sided, and test statistic is a simple log fold change
#' 
#' A paired permutation test is implemented specifically for this package. The test is similar to the original, but with a different test statistic and permutation scheme. The permutations are constrained in the paired version such that the predictor is only permuted within each level of the paired argument (e.g. subjects). The test statistic first finds the log-ratio between the two predictor levels (e.g. case and control) for each level of the paired argument and the final statistic is the mean of these log-ratios.
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param testStat Function. Function for the test statistic. Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param testStat.pair Function. Function for test statistc for paired analysis. Should take two vectors as arguments. Default is a log fold change: mean(log((case abundances+1)/(control abundances+1)))
#' @param noOfIterations Integer. Iterations for permutations. Default 10000
#' @param margin Numeric. Margin for when to stop iterations if p-value is high and unlikely to become low
#' @export

DA.per <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, noOfIterations = 10000, margin = 50){

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    } else {
    count_table <- data
  }
  
  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    predictor <- predictor[order(paired)]
    testStat <- testStat.pair
  }
  
  predictor <- as.numeric(as.factor(predictor))-1
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  nullStatList <- list()
  
  # Create shuffled predictors
  shuffledpredictorsList <- list()
  if(is.null(paired)){
    for (k in 1:noOfIterations){
      shuffledpredictorsList[[k]] <- sample(predictor)
    }
  } else {
    for (k in 1:noOfIterations){
      shuffledpredictorsList[[k]] <- unlist(lapply(1:(length(predictor)/2),function(x) sample(c(0,1))))
    }
  }
  
  iterations <- nrow(count_table)
  p <- numeric(iterations)
  FC <- numeric(iterations)
  coverage <- numeric(iterations)
  
  for(i in 1:iterations){
    
    count_row      <- as.numeric(count_table[i,])
    
    real_case    <- count_row[predictor==1]
    real_control <- count_row[predictor==0]
    realStat     <- testStat(real_case,real_control) 
    
    Wnull <- numeric(noOfIterations)
    
    for(j in 1:noOfIterations){
      case     <- count_row[shuffledpredictorsList[[j]]==1]
      control  <- count_row[shuffledpredictorsList[[j]]==0]
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
  output_df$Method <- "Permutation (per)"
  
  if(class(data) == "phyloseq") output_df <- add.tax.DA(data, output_df)
  
  return(output_df)
}


