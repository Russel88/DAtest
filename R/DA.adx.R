#' Aldex t.test and wilcox
#' 
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param ... Additional arguments for the aldex function
#' @export

DA.adx <- function(count_table, outcome, ...){
  
  library(ALDEx2, quietly = TRUE)
  
  x <- aldex(data.frame(count_table), outcome, verbose = FALSE, ...)
  x$Feature <- rownames(x)
  return(x)
  
}




