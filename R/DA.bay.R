#' baySeq
#' 
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the getLikelihoods function
#' @param ... Additional arguments to the getPriors.NB and getLikelihoods functions
#' @export

DA.bay <- function(data, predictor, p.adj = "fdr", allResults = FALSE, ...){
  
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
  
  DA.bay.args <- list(...)
  getPriors.NB.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(getPriors.NB))]
  getLikelihoods.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(getLikelihoods))]
  
  CD <- do.call(getPriors.NB, c(list(cD=CD, cl=NULL), getPriors.NB.args))
  CD <- do.call(getLikelihoods, c(list(cD=CD, cl=NULL), getLikelihoods.args))
  
  tc <- topCounts(CD, group = "DE", number=nrow(count_table))
  tc <- tc[,c(1,rev(ncol(tc)-0:4))]
  
  if(is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj))
  if(!is.null(tc$ordering))
    output_df <- data.frame(Feature = as.character(tc$annotation), pval = 1 - tc$Likelihood, pval.adj = p.adjust(1 - tc$Likelihood, method = p.adj), ordering = tc$ordering)
  
  
  output_df$Method <- "baySeq (bay)"

  if(class(data) == "phyloseq") output_df <- add.tax.DA(data, output_df)
  
  if(allResults) return(CD) else return(output_df)
}

