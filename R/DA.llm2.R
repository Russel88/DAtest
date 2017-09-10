#' Linear regression
#'
#' With log transformation of relative abundances.
#' Mixed-effect model is used when a paired argument is included, with the paired variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 0.001
#' @param ... Additional arguments for the lm/lme functions
#' @export

DA.llm2 <- function(data, predictor, paired = NULL, p.adj = "fdr", delta = 0.001, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    } else {
    count_table <- data
  }
  
  count_table <- as.data.frame.matrix(count_table)
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  if(is.null(paired)){
    lmr <- function(x){
      fit <- NULL
      tryCatch(
        fit <- lm(x ~ predictor, ...), 
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
        fit <- lme(x ~ predictor, random = ~1|paired, ...), 
        error = function(x) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,c(1,2,4,5)]
        } else NA
      } else NA 
    }
  }
  
  res <- as.data.frame(t(as.data.frame(apply(count.rel,1,lmr))))
  colnames(res) <- c("Estimate","Std.Error","t-value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Log Linear regression 2"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  return(res)
  
}
