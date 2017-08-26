#' RAIDA
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param predictor Factor. The outcome of interest. E.g. case and control
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the raida function
#' @export

DA.rai <- function(count_table, outcome, p.adj = "fdr", ...){
  
  library(RAIDA, quietly = TRUE)
  
  count_table.o <- as.data.frame(count_table[,order(outcome)])
  
  res <- raida(count_table.o, n.lib = as.numeric(table(outcome)), mtcm = p.adj, ...)
  res$Feature <- rownames(res)
  colnames(res)[1] <- "pval"
  colnames(res)[2] <- "pval.adj"
  res$Method <- "RAIDA" 
  
  return(res)

}




