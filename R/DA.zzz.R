#' User-defined function
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param dafun Function to apply to data. Should take input in the following order: count_table (data.frame, samples are columns), predictor (vector), paired (factor), covars (named list with vector).
#' @export
DA.zzz <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", dafun = NULL){
  
  if(is.null(dafun)) stop("dafun has to be defined")
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    if(!is.null(covars)){
      covars.n <- covars
      covars <- list()
      for(i in 1:length(covars.n)){
        covars[[i]] <- suppressWarnings(as.matrix(sample_data(data)[,covars.n[i]]))
      }
      names(covars) <- covars.n
    } 
  } else {
    count_table <- data
  }

  res <- dafun(count_table, predictor, paired, covars)
  
  # Check output
  if(!all(c("pval","Feature","Method") %in% colnames(res))) stop("The following columns has to be present in output: pval, Feature and Method")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  
  return(res)
}




