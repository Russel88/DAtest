#' Extract data from a \code{phyloseq} object to be used in \code{DAtest}
#' 
#' @param data A \code{phyloseq} object
#' @param predictor The \code{predictor} of interest
#' @param paired Factor for paired/blocked experimental designs.
#' @param covars A character vector with names of the variables in \code{sample_data(data)}
#' @export

DA.phyloseq <- function(data, predictor, paired = NULL, covars = NULL){
 
  if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
  if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
  if(!is.null(paired)){
    if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
  }
  count_table <- otu_table(data)
  if(!taxa_are_rows(data)) count_table <- t(count_table)
  predictor <- unlist(sample_data(data)[,predictor])
  if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
  if(!is.null(covars)){
    covars.n <- covars
    covars <- list()
    for(i in 1:length(covars.n)){
      covars[[i]] <- unlist(sample_data(data)[,covars.n[i]])
    }
    names(covars) <- covars.n
  } 

  DA.data <- list(count_table = count_table,
                  predictor = predictor,
                  paired = paired,
                  covars = covars)
  
  return(DA.data)
  
}
