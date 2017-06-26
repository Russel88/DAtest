#' ANCOM
#'
#' @export

DA.anc <- function(otu_table, outcome, sig = 0.05, multcorr = 3, tau = 0.02, theta = 0.1, repeated = FALSE){
  
  library(ancom.R, quietly = TRUE)

  input <- as.data.frame(t(otu_table))
  input$Group <- outcome
  
  res <- ANCOM(input, sig, multcorr, tau, theta, repeated) 

  df <- data.frame(OTU = rownames(otu_table),
                   W = res[[1]],
                   pval = 1,
                   pval.adj = 1)  
  df[df$OTU %in% res[[2]],"pval"] <- 0
  df[df$OTU %in% res[[2]],"pval.adj"] <- 0
  
  df$Method <- "ANCOM"
  
  return(df)  
}
