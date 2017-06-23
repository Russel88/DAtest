#' Spike-in
#'
#' Modified version of the one from:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

spikein <- function(otu_table, outcome, spikeMethod = "mult", effectSize = 2, k){
  
  if(effectSize == 1) spikeMethod <- "none"

  otu_table <- as.data.frame(otu_table)
  outcome <- as.numeric(as.factor(outcome))-1
  
  # Choose OTUs to spike
  propotu <- apply(otu_table,2,function(x) x/sum(x))
  otu_abundances <- sort(rowSums(propotu)/ncol(propotu))

  # Only spike OTUs present in cases
  approved_otu_abundances <- otu_abundances[ 
    names(otu_abundances) %in% row.names( otu_table[ rowSums(otu_table[,outcome == 1]) > 0, outcome == 1] ) ]
  lower_tert <- names(approved_otu_abundances[approved_otu_abundances < quantile(approved_otu_abundances,1/3)])
  mid_tert <- names(approved_otu_abundances[approved_otu_abundances >= quantile(approved_otu_abundances,1/3) & approved_otu_abundances < quantile(approved_otu_abundances,2/3)])
  upper_tert <- names(approved_otu_abundances[approved_otu_abundances >= quantile(approved_otu_abundances,2/3)])
  spiked_otus <- c( 	sample(lower_tert, k)	,
                     sample(mid_tert, k)		,
                     sample(upper_tert,k)	)
  spiked_otu_index <- which(row.names(otu_table) %in% spiked_otus)
  
  # Spike OTUs either by addition or multiplication
  oldSums <- colSums(otu_table)

  if(spikeMethod == "mult"){
    otu_table[spiked_otu_index,outcome==1] <- otu_table[spiked_otu_index, outcome==1] * effectSize
  }
  
  if(spikeMethod == "add"){
    nonzeroMeans <- lapply(spiked_otu_index, function(j){
      v <- propotu[j,]
      v <- v[v != 0]
      mean(v)
    })
    cases <- otu_table[,outcome == 1]
    caseDepth <- colSums(cases)
    addlist <- lapply(nonzeroMeans, function(l) l*caseDepth*effectSize)
    addcounts <- do.call(rbind, addlist)
    addcounts[cases[spiked_otu_index,] == 0] <- 0
    cases[spiked_otu_index,] <- cases[spiked_otu_index,]+addcounts
    otu_table[,outcome == 1] <- cases	
  }
  
  newSums <- colSums(otu_table)
  otu_table<- round(t(t(otu_table) * oldSums/newSums))
  list(otu_table,spiked_otus)
  
}