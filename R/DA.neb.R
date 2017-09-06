#' Negative binomial glm
#'
#' With log(librarySize) as offset.
#' Mixed-effect model is used when a paired argument is included, with the paired variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param outcome The outcome of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the glm.nb/glmer.nb functions
#' @import MASS
#' @importFrom lme4 glmer.nb glmer
#' @export

DA.neb <- function(data, outcome, paired = NULL, p.adj = "fdr", ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(outcome) > 1 | length(paired) > 1) stop("When data is a phyloseq object outcome and paired should only contain the name of the variables in sample_data")
    if(!outcome %in% sample_variables(data)) stop(paste(outcome,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    outcome <- suppressWarnings(as.matrix(sample_data(data)[,outcome]))
    if(!is.null(paired)) paired <- suppressWarnings(as.matrix(sample_data(data)[,paired]))
  } else {
    count_table <- data
  }
  
  libSize <- colSums(count_table)
  count_table <- as.data.frame.matrix(count_table)
  
  if(is.null(paired)){
    negbin <- function(x){
      fit <- NULL
      tryCatch(
        fit <- glm.nb(x ~ outcome + offset(log(libSize)), ...), 
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
        fit <- lme4::glmer.nb(x ~ outcome + offset(log(libSize)) + (1|paired), ...), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,]
        } else NA
      } else NA 
    }
  }
  
  res <- as.data.frame(t(as.data.frame(apply(count_table,1,negbin))))
  colnames(res) <- c("Estimate","Std.Error","z value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Negbinom"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  return(res)
  
}
