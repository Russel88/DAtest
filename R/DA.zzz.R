#' User-defined function
#' 
#' Apply a user-defined function to all features of a count table. For implemetation in \code{testDA} and \code{allDA}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param FUN Function to apply to data. Should take input in the following order: \code{count_table} (data.frame, samples are columns), \code{predictor} (vector), \code{paired} (factor), \code{covars} (named list with vector). Output should be a dataframe with at least the following columns: \code{Feature}, \code{pval} and \code{Method}.
#' @return A data.frame with results from the user-defined method
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 100, ncol = 10)
#' rownames(mat) <- 1:100
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Define function for t-test
#' myfun <- function(count_table, predictor, paired, covars){ 
#'  
#' # Relative abundance
#' rel <- apply(count_table, 2, function(x) x/sum(x))
#' 
#' # t-test function
#' # Wrapping this function in tryCatch(..., error = function(e){NA}) 
#' # ensures that our main function won't fail if t.test fails on some features
#' tfun <- function(x){
#'     tryCatch(t.test(x ~ predictor)$p.value, error = function(e){NA}) 
#' }
#' 
#' # P-values for each feature
#' pvals <- apply(rel, 1, tfun)
#' 
#' # Collect and return data
#' df <- data.frame(Feature = rownames(count_table),
#'                  pval = pvals)
#' df$pval.adj <- p.adjust(df$pval, method = "fdr")
#' df$Method <- "My own t-test"
#' return(df)
#' }
#' 
#' # Running the test
#' res <- DA.zzz(data = mat, predictor = pred, FUN = myfun)
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
    for(i in seq_along(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }

  # Run the test
  res <- FUN(count_table, predictor, paired, covars)
  
  # Check output
  if(!all(c("pval","Feature","Method") %in% colnames(res))) stop("The following columns has to be present in output: pval, Feature and Method")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
  if(class(data) == "phyloseq") res <- addTax(data, res)
  
  return(res)
}




