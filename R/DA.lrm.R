#' Linear regression
#' 
#' Mixed-effect model is a paired argument is present, with the paired variable as random intercept
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the lm/lme functions
#' @import nlme
#' @export

DA.lrm <- function(count_table, outcome, paired = NULL, relative = TRUE, p.adj = "fdr", ...){
 
  library(nlme)
  
  count_table <- as.data.frame.matrix(count_table)
  if(relative) count_table <- apply(count_table,2,function(x) x/sum(x))
  
  if(is.null(paired)){
    lmr <- function(x){
      fit <- NULL
      tryCatch(
        fit <- lm(x ~ outcome, ...), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,]
        } else NA
      } else NA 
    }
  } else {
    lmr <- function(x){
      fit <- NULL
      tryCatch(
        fit <- lme(x ~ outcome, random = ~1|paired, ...), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,c(1,2,4,5)]
        } else NA
      } else NA 
    }
  }
  
  res <- as.data.frame(t(as.data.frame(apply(count_table,1,lmr))))
  colnames(res) <- c("Estimate","Std.Error","t-value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Linear regression"
  return(res)
  
}
