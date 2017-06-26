#' baySeq
#'
#' From
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.bay <- function(otu_table, outcome){
  
  library(baySeq, quietly = TRUE)
  
  outcome <- as.numeric(as.factor(outcome))-1
  CD <- new("countData", data=as.matrix(otu_table), replicates = ifelse(as.logical(outcome), "simA", "simB"), groups = list(NDE = rep(1,length(outcome)),DE=ifelse(as.logical(outcome),1,2))) # simA = cases
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- data.frame(name=rownames(otu_table))
  
  CD <- getPriors.NB(CD, cl= NULL)
  CD <- getLikelihoods(CD, cl = NULL)
  
  tc <- topCounts(CD, group = "DE", number=nrow(otu_table))
  tc <- tc[,c(1,rev(ncol(tc)-0:4))]
  
  if(is.null(tc$ordering))
    output_df <- data.frame(OTU = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = tc$FDR.DE)
  if(!is.null(tc$ordering))
    output_df <- data.frame(OTU = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = tc$FDR.DE, ordering = tc$ordering)
  
  
  output_df$Method <- "baySeq"

  return(output_df)
}

