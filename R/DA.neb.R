#' Negative binomial glm
#'
#' @import MASS
#' @export

DA.neb <- function(otu_table, outcome, p.adj){
 
  library(MASS, quietly = TRUE)
  
  libSize <- colSums(otu_table)
  otu_table <- as.data.frame.matrix(otu_table)
  
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
  res <- as.data.frame(t(as.data.frame(apply(otu_table,1,negbin))))
  colnames(res) <- c("Estimate","Std.Error","z value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$OTU <- rownames(res)
  res$Method <- "Negbinom"
  return(res)
  
}
