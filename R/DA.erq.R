#' EdgeR quasi-likelihood
#' 
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the glmQLFit function
#' @export

DA.erq <- function(count_table, outcome, paired = NULL, p.adj = "fdr", ...){
  
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
  fit <- glmQLFit(y,design, ...)
  if(is.numeric(outcome)){
    qlf <- glmQLFTest(fit,coef=2)
    ta <- qlf$table
    colnames(ta)[4] <- "pval"
  } else {
    qlf <- glmQLFTest(fit,coef=seq(2,length(levels(as.factor(outcome)))))
    ta <- qlf$table
    colnames(ta)[(length(levels(as.factor(outcome)))+2)] <- "pval"
  }

  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR qll"
  
  return(ta)
  
}


