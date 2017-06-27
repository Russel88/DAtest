#' Negative binomial glm
#'
#' @import MASS lme4
#' @export

DA.neb <- function(otu_table, outcome, paired = NULL, p.adj){
 
  library(MASS, quietly = TRUE)
  library(lme4, quietly = TRUE)
  
  libSize <- colSums(otu_table)
  otu_table <- as.data.frame.matrix(otu_table)
  
  if(is.null(paired)){
    negbin <- function(x){
      fit <- NULL
      tryCatch(
        fit <- glm.nb(x ~ outcome + offset(log(libSize))), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,]
        } else NA
      } else NA 
    }
  } else {
    negbin <- function(x){
      fit <- NULL
      tryCatch(
        fit <- glmer.nb(x ~ outcome + offset(log(libSize)) + (1|paired)), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,]
        } else NA
      } else NA 
    }
  }
  
  res <- as.data.frame(t(as.data.frame(apply(otu_table,1,negbin))))
  colnames(res) <- c("Estimate","Std.Error","z value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$OTU <- rownames(res)
  res$Method <- "Negbinom"
  return(res)
  
}
