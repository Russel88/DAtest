#' baySeq
#' 
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param samplesize How large a sample should be taken in estimating the priors? Default 1e5
#' @param samplingSubset If given, the priors will be sampled only from the subset specified. Default NULL
#' @param equalDispersions Should we assume equal dispersions of data across all groups in the 'cD' object? Defaults to TRUE
#' @param estimation Defaults to "QL", indicating quasi-likelihood estimation of priors. Currently, the only other possibilities are "ML", a maximum-likelihood method, and "edgeR", the moderated dispersion estimates produced by the 'edgeR' package
#' @param zeroML Should parameters from zero data (rows that within a group are all zeros) be estimated using maximum likelihood methods (which will result in zeros in the parameters?
#' @param consensus If TRUE, creates a consensus distribution rather than a separate distribution for each member of the groups structure in the 'cD' object
#' @param ... Additional arguments to the getLikelihoods function
#' @export

DA.bay <- function(count_table, outcome, p.adj = "fdr", samplesize = 1e5, samplingSubset = NULL, equalDispersions = TRUE, estimation = "QL", zeroML = FALSE, consensus = FALSE, ...){
  
  library(baySeq, quietly = TRUE)
  
  outcome <- as.numeric(as.factor(outcome))-1
  CD <- new("countData", data=as.matrix(count_table), replicates = ifelse(as.logical(outcome), "simA", "simB"), groups = list(NDE = rep(1,length(outcome)),DE=ifelse(as.logical(outcome),1,2))) # simA = cases
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- data.frame(name=rownames(count_table))
  
  CD <- getPriors.NB(CD, cl= NULL, samplesize=samplesize, samplingSubset=samplingSubset, equalDispersions=equalDispersions, estimation=estimation, zeroML=zeroML, consensus=consensus)
  CD <- getLikelihoods(CD, cl = NULL, ...)
  
  tc <- topCounts(CD, group = "DE", number=nrow(count_table))
  tc <- tc[,c(1,rev(ncol(tc)-0:4))]
  
  if(is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj))
  if(!is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj), ordering = tc$ordering)
  
  
  output_df$Method <- "baySeq"

  return(output_df)
}

