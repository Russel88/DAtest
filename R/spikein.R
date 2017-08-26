#' Spike-in
#'
#' Internal function for the testDA function.
#' 
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' 
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor or Numeric. The outcome of interest. E.g. case and control. If the predictor has more than two levels, only the 2. level will be spiked. If the predictor is numeric it will be treated as such in the analyses
#' @param spikeMethod Character. Multiplicative ("mult") or additive ("add") spike-in. Default "mult". Use "add" if there are negative values in count_table
#' @param effectSize Integer. The effect size for the spike-ins. Default 2
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). k=c(5,10,15): 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default c(5,5,5)
#' @param num.pred Logical. Is the outcome numeric? Default FALSE
#' @export

spikein <- function(count_table, outcome, spikeMethod = "mult", effectSize = 2, k, num.pred = FALSE){
  
  if(num.pred == TRUE & spikeMethod == "add") stop("Cannot use additive spike-in if predictor is numeric")
  if(effectSize < 0) stop("Effect size should be positive")
  if(effectSize == 1) spikeMethod <- "none"
  
  count_table <- as.data.frame(count_table)
  outcome <- as.numeric(as.factor(outcome))-1
  
  # Choose Features to spike
  propcount <- apply(count_table,2,function(x) x/sum(x))
  count_abundances <- sort(rowSums(propcount)/ncol(propcount))

  # Only spike Features present in cases
  approved_count_abundances <- count_abundances[ 
    names(count_abundances) %in% row.names( count_table[ rowSums(count_table[,outcome == 1]) > 0, outcome == 1] ) ]
  
  lower_tert <- names(approved_count_abundances[approved_count_abundances < quantile(approved_count_abundances,1/3)])
  mid_tert <- names(approved_count_abundances[approved_count_abundances >= quantile(approved_count_abundances,1/3) & approved_count_abundances < quantile(approved_count_abundances,2/3)])
  upper_tert <- names(approved_count_abundances[approved_count_abundances >= quantile(approved_count_abundances,2/3)])
  
  spike_features <- c(sample(lower_tert, k[1]), sample(mid_tert, k[2]), sample(upper_tert,k[3]))
  spike_feature_index <- which(row.names(count_table) %in% spike_features)
  
  # Spike Features either by addition or multiplication
  oldSums <- colSums(count_table)

  if(spikeMethod == "mult"){
    if(num.pred){
      outcome <- as.numeric(outcome)
      count_table[spike_feature_index,] <- count_table[spike_feature_index, ] * (effectSize ^ outcome)
    } else {
      count_table[spike_feature_index,outcome==1] <- count_table[spike_feature_index, outcome==1] * effectSize
    }
  }
  
  if(spikeMethod == "add"){
    nonzeroMeans <- lapply(spike_feature_index, function(j){
      v <- propcount[j,]
      v <- v[v != 0]
      mean(v)
    })
    cases <- count_table[,outcome == 1]
    caseDepth <- colSums(cases)
    addlist <- lapply(nonzeroMeans, function(l) l*caseDepth*effectSize)
    addcounts <- do.call(rbind, addlist)
    addcounts[cases[spike_feature_index,] == 0] <- 0
    cases[spike_feature_index,] <- cases[spike_feature_index,]+addcounts
    count_table[,outcome == 1] <- cases	
  }
  
  newSums <- colSums(count_table)
  count_table <- round(t(t(count_table) * oldSums/newSums))
  list(count_table,spike_features)
  
}