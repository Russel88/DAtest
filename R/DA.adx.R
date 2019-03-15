#' Aldex t.test and wilcox
#' 
#' Implementation of \code{aldex} for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param ... Additional arguments for the \code{aldex} function
#' @return A data.frame with with results.
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 100, ncol = 10)
#' rownames(mat) <- 1:100
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running ALDEx2
#' res <- DA.adx(data = mat, predictor = pred)
#' @export

DA.adx <- function(data, predictor, ...){
  
  ok <- tryCatch({
    loadNamespace("ALDEx2")
    TRUE
  }, error=function(...) FALSE)
  
  if (ok) {

    # Extract from phyloseq
    if(class(data) == "phyloseq"){
      DAdata <- DA.phyloseq(data, predictor)
      count_table <- DAdata$count_table
      predictor <- DAdata$predictor
    } else {
      count_table <- data
    }
    
    # Run test
    x <- ALDEx2::aldex(data.frame(count_table), as.character(predictor), ...)
    x$ordering <- NA
    x[!is.na(x$effect) & x$effect > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
    x[!is.na(x$effect) & x$effect < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
    x$Feature <- rownames(x)
    
    if(class(data) == "phyloseq") x <- addTax(data, x)
    
    return(x)
  } else {
    stop("ALDEx2 package required")
  }
  
}




