#' Aldex t.test and wilcox

#' @export

DA.adx <- function(count_table, outcome, mc.samples = 128, p.adj){
  
  library(ALDEx2, quietly = TRUE)
  
  x <- aldex(data.frame(count_table), outcome, mc.samples = mc.samples, verbose = FALSE)
  x$Feature <- rownames(x)
  return(x)
  
}




