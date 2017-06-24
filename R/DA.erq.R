#' EdgeR quasi-likelihood

#' @export

DA.erq <- function(otu_table, outcome, p.adj){
  
  library(edgeR, quietly = TRUE)
  otu_table <- as.data.frame(otu_table)
  y <- DGEList(counts=otu_table,genes = data.frame(OTU = row.names(otu_table)))
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~outcome)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  ta <- qlf$table
  colnames(ta)[4] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$OTU <- rownames(ta)
  ta$Method <- "EdgeR qll"
  
  return(ta)
  
}


