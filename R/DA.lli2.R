#' LIMMA
#' 
#' With log transformation of relative abundances.
#' Some implementation is borrowed from:
#' http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for log transformation. Default 0.001
#' @export

DA.lli2 <- function(count_table, outcome, paired = NULL, p.adj = "fdr", delta = 0.001, ...){
  
  library(limma, quietly = TRUE)
  
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  
  count.rel <- log(count.rel + delta)
    
  if(is.null(paired)) design <- model.matrix(~outcome) else design <- model.matrix(~as.factor(paired)+outcome)
  n <- dim(count.rel)[1]
  fit <- lmFit(count.rel, design)
  fit.eb <- eBayes(fit, ...)
  if(is.null(paired)) Estimate <- fit.eb$coefficients else Estimate <- fit.eb$coefficients[,c(1,(length(levels(as.factor(paired)))+1):ncol(fit.eb$coefficients))]
  df.residual <- fit.eb$df.residual
  df.prior <- rep(fit.eb$df.prior, n)
  s2.prior <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.stat <- fit.eb$t[, 2]
  pval <- fit.eb$p.value[, 2]
  pval.adj <- p.adjust(pval, method = p.adj)
  res <- data.frame(Estimate, t.stat, pval, pval.adj, df.residual, df.prior, s2.prior, s2, s2.post)
  res$Feature <- rownames(res)
  res$Method <- "Log LIMMA 2"

  return(res)  
}

