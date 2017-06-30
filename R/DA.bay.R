#' baySeq
#'
#' From
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.bay <- function(count_table, outcome, p.adj){
  
  library(baySeq, quietly = TRUE)
  
  outcome <- as.numeric(as.factor(outcome))-1
  CD <- new("countData", data=as.matrix(count_table), replicates = ifelse(as.logical(outcome), "simA", "simB"), groups = list(NDE = rep(1,length(outcome)),DE=ifelse(as.logical(outcome),1,2))) # simA = cases
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- data.frame(name=rownames(count_table))
  
  CD <- getPriors.NB(CD, cl= NULL)
  CD <- getLikelihoods(CD, cl = NULL)
  
  tc <- topCounts(CD, group = "DE", number=nrow(count_table))
  tc <- tc[,c(1,rev(ncol(tc)-0:4))]
  
  if(is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj))
  if(!is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj), ordering = tc$ordering)
  
  
  output_df$Method <- "baySeq"

  return(output_df)
}

