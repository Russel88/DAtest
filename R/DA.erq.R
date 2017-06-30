#' EdgeR quasi-likelihood

#' @export

DA.erq <- function(count_table, outcome, paired = NULL, p.adj){
  
  library(edgeR, quietly = TRUE)
  otu_table <- as.data.frame(count_table)
  y <- DGEList(counts=count_table,genes = data.frame(Feature = row.names(count_table)))
  y <- edgeR::calcNormFactors(y)
  if(is.null(paired)){
    design <- model.matrix(~outcome)
  } else {
    design <- model.matrix(~outcome + paired)
  }
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  ta <- qlf$table
  colnames(ta)[4] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR qll"
  
  return(ta)
  
}


