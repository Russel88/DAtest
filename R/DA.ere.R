#' EdgeR
#'
#' From:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.ere <- function(otu_table, outcome){
  
  library(edgeR, quietly = TRUE)
  otu_table <- as.data.frame(otu_table)
  x <- DGEList(counts = otu_table, group = outcome, genes = data.frame(OTU = row.names(otu_table)))
  x <- edgeR::calcNormFactors(x)
  x <- estimateTagwiseDisp(estimateCommonDisp(x))
  ta <- exactTest(x)[[1]]
  colnames(ta)[3] <- "pval"
  ta$OTU <- rownames(ta)
  ta$Method <- "EdgeR exact"
 
  return(ta) 
}


