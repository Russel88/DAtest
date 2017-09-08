#' baySeq
#' 
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param samplesize How large a sample should be taken in estimating the priors? Default 1e5
#' @param samplingSubset If given, the priors will be sampled only from the subset specified. Default NULL
#' @param equalDispersions Should we assume equal dispersions of data across all groups in the 'cD' object? Defaults to TRUE
#' @param estimation Defaults to "QL", indicating quasi-likelihood estimation of priors. Currently, the only other possibilities are "ML", a maximum-likelihood method, and "edgeR", the moderated dispersion estimates produced by the 'edgeR' package
#' @param zeroML Should parameters from zero data (rows that within a group are all zeros) be estimated using maximum likelihood methods (which will result in zeros in the parameters?
#' @param consensus If TRUE, creates a consensus distribution rather than a separate distribution for each member of the groups structure in the 'cD' object
#' @param ... Additional arguments to the getLikelihoods function
#' @export

DA.bay <- function(data, predictor, p.adj = "fdr", samplesize = 1e5, samplingSubset = NULL, equalDispersions = TRUE, estimation = "QL", zeroML = FALSE, consensus = FALSE, ...){
  
  library(baySeq)
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
  } else {
    count_table <- data
  }
  
  predictor <- as.numeric(as.factor(predictor))-1
  CD <- new("countData", data=as.matrix(count_table), replicates = ifelse(as.logical(predictor), "simA", "simB"), groups = list(NDE = rep(1,length(predictor)),DE=ifelse(as.logical(predictor),1,2))) # simA = cases
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

  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      output_df <- merge(output_df, tax, by.x = "Feature", by.y = "row.names")
      rownames(output_df) <- NULL
    } 
  }
  
  return(output_df)
}

