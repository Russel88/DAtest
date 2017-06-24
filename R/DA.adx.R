#' Aldex t.test and wilcox

#' @export

DA.adx <- function(otu_table, outcome, mc.samples = 128, p.adj){
  
  library(ALDEx2, quietly = TRUE)
  
  x <- aldex(data.frame(otu_table), outcome, mc.samples = mc.samples, verbose = FALSE)
  x$OTU <- rownames(x)
  return(x)
  
}




