#' ANCOM
#'
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param ... Additional arguments for the ANCOM function
#' @export

DA.anc <- function(count_table, outcome, ...){
  
  library(ancom.R, quietly = TRUE)

  input <- as.data.frame(t(count_table))
  input$Group <- outcome
  
  res <- ANCOM(input, ...) 

  df <- data.frame(Feature = rownames(count_table),
                   W = res[[1]],
                   pval = 1,
                   pval.adj = 1)  
  df[df$Feature %in% res[[2]],"pval"] <- 0
  df[df$Feature %in% res[[2]],"pval.adj"] <- 0
  
  df$Method <- "ANCOM"
  
  return(df)  
}
