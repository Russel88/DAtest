#' baySeq
#' 
#' Implementation of baySeq for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param allResults If TRUE will return raw results from the \code{getLikelihoods} function
#' @param ... Additional arguments to the \code{getPriors.NB} and \code{getLikelihoods} functions
#' @return A data.frame with with results.
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(100, size = 0.1, mu = 500), nrow = 50, ncol = 10)
#' rownames(mat) <- 1:50
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running baySeq
#' res <- DA.bay(data = mat, predictor = pred)
#' @export

DA.bay <- function(data, predictor, allResults = FALSE, ...){
  
  ok <- tryCatch({
    loadNamespace("baySeq")
    TRUE
  }, error=function(...) FALSE)
  
  if (ok){
    
    # Extract from phyloseq
    if(is(data, "phyloseq")){
      DAdata <- DA.phyloseq(data, predictor)
      count_table <- DAdata$count_table
      predictor <- DAdata$predictor
    } else {
      count_table <- data
    }
    
    # Collect data
    CD <- new("countData", data=as.matrix(count_table), replicates = ifelse(as.logical(as.numeric(as.factor(predictor))-1), "simA", "simB"), groups = list(NDE = rep(1,length(predictor)),DE=predictor)) # simA = cases
    baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
    CD@annotation <- data.frame(name=rownames(count_table))
    
    # Arguments
    DA.bay.args <- list(...)
    getPriors.NB.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(baySeq::getPriors.NB))]
    getLikelihoods.args <- DA.bay.args[names(DA.bay.args) %in% names(formals(baySeq::getLikelihoods))]
    
    # Run test
    CD <- do.call(baySeq::getPriors.NB, c(list(cD=CD, cl=NULL), getPriors.NB.args))
    CD <- do.call(baySeq::getLikelihoods, c(list(cD=CD, cl=NULL), getLikelihoods.args))
    
    # Extract results
    tc <- baySeq::topCounts(CD, group = "DE", number=nrow(count_table))
    tc <- tc[,c(1,rev(ncol(tc)-0:3))]
    
    output_df <- data.frame(Feature = as.character(tc$name), pval = (1 - tc$likes), pval.adj = tc$FDR.DE, ordering = tc$DE)
    
    output_df$Method <- "baySeq (bay)"
    
    if(is(data, "phyloseq")) output_df <- addTax(data, output_df)
    
    if(allResults) return(CD) else return(output_df)
  } else {
    stop("baySeq package required")
  }
  
}

