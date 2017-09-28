#' Linear regression
#' 
#' Mixed-effect model is used when a paired argument is present, with the paired variable as random intercept
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param out.anova If TRUE will output results and p-values from anova. If false will output results for 2. level of the predictor.
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the lm/lme function
#' @param ... Additional arguments for the lm/lme functions
#' @import nlme
#' @export

DA.lrm <- function(data, predictor, paired = NULL, covars = NULL, relative = TRUE, out.anova = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
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
    if(!is.null(covars)){
      for(i in 1:length(covars)){
        assign(covars[i], suppressWarnings(as.matrix(sample_data(data)[,covars[i]])))
      }
    } 
  } else {
    count_table <- data
    if(!is.null(covars)){
      for(i in 1:length(covars)){
        assign(names(covars)[i], covars[[i]])
      }
    }
  }
  
  count_table <- as.data.frame.matrix(count_table)
  if(relative) count_table <- apply(count_table,2,function(x) x/sum(x))
  
  if(is.null(covars)){
    form <- paste("x ~ predictor")
  } else {
    if(class(data) == "phyloseq"){
      form <- paste("x ~ predictor+",paste(covars, collapse="+"),sep = "")
    } else {
      form <- paste("x ~ predictor+",paste(names(covars), collapse="+"),sep = "")
    }
  }
  
  if(is.null(paired)){
    lmr <- function(x){
      fit <- NULL
      tryCatch(
        fit <- lm(as.formula(form), ...), 
        error = function(e) fit <- NULL)
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
        fit <- lme(as.formula(form), random = ~1|paired, ...), 
        error = function(e) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          coef(summary(fit))[2,c(1,2,4,5)]
        } else NA
      } else NA 
    }
  }
  
  if(out.anova){
    if(is.null(paired)){
      lmr <- function(x){
        fit <- NULL
        tryCatch(
          fit <- anova(lm(as.formula(form),...))[1,], 
          error = function(e) fit <- NULL)
      }
    } else {
      lmr <- function(x){
        fit <- NULL
        tryCatch(
          fit <- anova(lme(as.formula(form), random = ~1|paired,...))[2,], 
          error = function(e) fit <- NULL)
      }
    }
  }
  
  if(out.anova){
    if(is.null(paired)){
      res <- as.data.frame(do.call(rbind,apply(count_table,1,lmr)))
      colnames(res) <- c("Df","Sum Sq","Mean Sq","F value","pval")
    } else {
      res <- as.data.frame(do.call(rbind,apply(count_table,1,lmr)))
      colnames(res) <- c("numDF","denDF","F-value","pval")
    }
  } else {
    res <- as.data.frame(t(as.data.frame(apply(count_table,1,lmr))))
    colnames(res) <- c("Estimate","Std.Error","t-value","pval")
  }
  
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Linear regression (lrm)"
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  
  if(allResults){
    if(is.null(paired)){
      lmr <- function(x){
        fit <- NULL
        tryCatch(fit <- lm(as.formula(form), ...), error = function(e) fit <- NULL)  
      }
    } else {
      lmr <- function(x){
        fit <- NULL
        tryCatch(
          fit <- lme(as.formula(form), random = ~1|paired, ...), error = function(e) fit <- NULL)
      }
    }
    return(apply(count_table,1,lmr))
  } else {
    return(res)
  }
  
}
