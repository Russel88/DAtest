#' EdgeR
#'
#' From:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.ere <- function(count_table, outcome, p.adj){
  
  library(edgeR, quietly = TRUE)
  otu_table <- as.data.frame(count_table)
  x <- DGEList(counts = count_table, group = outcome, genes = data.frame(Feature = row.names(count_table)))
  x <- edgeR::calcNormFactors(x)
  x <- estimateTagwiseDisp(estimateCommonDisp(x))
  ta <- exactTest(x)[[1]]
  colnames(ta)[3] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR exact"
 
  return(ta) 
}


