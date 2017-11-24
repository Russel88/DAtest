#' User-defined function
#' 
#' Apply a user-defined function to all features of a count table. For implemetation in \code{testDA} and \code{allDA}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param FUN Function to apply to data. Should take input in the following order: \code{count_table} (data.frame, samples are columns), \code{predictor} (vector), \code{paired} (factor), \code{covars} (named list with vector). Output should be a dataframe with at least the following columns: \code{Feature}, \code{pval} and \code{Method}.
#' @export
DA.zzz <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", FUN = NULL){
  
  if(is.null(FUN)) stop("FUN has to be defined")
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
    covars <- DAdata$covars
  } else {
    count_table <- data
  }
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }

  # Run the test
  res <- FUN(count_table, predictor, paired, covars)
  
  # Check output
  if(!all(c("pval","Feature","Method") %in% colnames(res))) stop("The following columns has to be present in output: pval, Feature and Method")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  
  return(res)
}




