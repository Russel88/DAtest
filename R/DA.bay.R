#' baySeq
#' 
#' Implementation of baySeq for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param allResults If TRUE will return raw results from the \code{getLikelihoods} function
#' @param ... Additional arguments to the \code{getPriors.NB} and \code{getLikelihoods} functions
#' @export

DA.bay <- function(data, predictor, p.adj = "fdr", allResults = FALSE, ...){
  
  suppressMessages(library(baySeq))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
  } else {
    count_table <- data
  }
  
  # Collect data
  CD <- new("countData", data=as.matrix(count_table), replicates = ifelse(as.logical(as.numeric(as.factor(predictor))-1), "simA", "simB"), groups = list(NDE = rep(1,length(predictor)),DE=predictor)) # simA = cases
  libsizes(CD) <- getLibsizes(CD)
  CD@annotation <- data.frame(name=rownames(count_table))
  
  # Arguments
  DA.bay.args <- list(...)
  getPriors.NB.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(getPriors.NB))]
  getLikelihoods.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(getLikelihoods))]
  
  # Run test
  CD <- do.call(getPriors.NB, c(list(cD=CD, cl=NULL), getPriors.NB.args))
  CD <- do.call(getLikelihoods, c(list(cD=CD, cl=NULL), getLikelihoods.args))
  
  # Extract results
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

