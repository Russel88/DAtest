#' Log linear regression
#'
#' Apply linear regression on multiple features with one \code{predictor}, with log transformation of relative abundances.
#' Mixed-effect model is used when a \code{paired} argument is included, with the \code{paired} variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in sample_data(data)
#' @param out.all If TRUE will output results and p-values from \code{anova}. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 0.001
#' @param coeff Integer. The p-value and log2FoldChange will be associated with this coefficient. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{lm}/\code{lme} function
#' @param ... Additional arguments for the \code{lm}/\code{lme} functions
#' @import nlme
#' @export

DA.llm2 <- function(data, predictor, paired = NULL, covars = NULL, out.all = NULL, p.adj = "fdr", delta = 0.001, coeff = 2, allResults = FALSE, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
    covars <- DAdata$covars
  } else {
    count_table <- data
  }
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }
  
  # out.all
  if(is.null(out.all)){
    if(length(unique(predictor)) == 2) out.all <- FALSE
    if(length(unique(predictor)) > 2) out.all <- TRUE
    if(is.numeric(predictor)) out.all <- FALSE
  }
  
  # Relative abundance and log
  count_table <- as.data.frame.matrix(count_table)
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  # Define design
  if(is.null(covars)){
    form <- paste("x ~ predictor")
  } else {
    form <- paste("x ~ predictor+",paste(names(covars), collapse="+"),sep = "")
  }
  
  # Define function
  if(is.null(paired)){
    lmr <- function(x){
      fit <- NULL
      tryCatch(
        fit <- lm(as.formula(form), ...), 
        error = function(e) fit <- NULL)
      if(!is.null(fit)) {
        if(nrow(coef(summary(fit))) > 1) {
          pval <- coef(summary(fit))[coeff,4]
          ests <- coef(summary(fit))[,1]
          c(ests,pval)
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
          pval <- coef(summary(fit))[coeff,5]
          ests <- coef(summary(fit))[,1]
          c(ests,pval)
        } else NA
      } else NA 
    }
  }
  
  ## for out.all TRUE
  if(out.all){
    if(is.null(paired)){
      lmr <- function(x){
        fit <- NULL
        tryCatch(
          fit <- lm(as.formula(form),...), 
          error = function(e) fit <- NULL)
        if(!is.null(fit)){
          ests <- coef(summary(fit))[,1]
          ano <- anova(fit)[1,]
          c(ano,ests)
        }
      }
    } else {
      lmr <- function(x){
        fit <- NULL
        tryCatch(
          fit <- lme(as.formula(form), random = ~1|paired, ...), 
          error = function(e) fit <- NULL)
        if(!is.null(fit)){
          ests <- coef(summary(fit))[,1]
          ano <- anova(fit)[2,]
          c(ano,ests)
        }
      }
    }
  }
  
  # Run tests
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
    if(out.all){
      if(is.null(paired)){
        res <- as.data.frame(do.call(rbind,apply(count_table,1,lmr)))
        colnames(res)[1:5] <- c("Df","Sum Sq","Mean Sq","F value","pval")
      } else {
        res <- as.data.frame(do.call(rbind,apply(count_table,1,lmr)))
        colnames(res)[1:4] <- c("numDF","denDF","F-value","pval")
      }
      res <- as.data.frame(lapply(res, unlist))
    } else {
      res <- as.data.frame(t(as.data.frame(apply(count_table,1,lmr))))
      colnames(res)[ncol(res)] <- "pval"
      res$log2FC <- log2((res[,1]+res[,coeff]) / res[,1])
      res[res[,1] < 0 & !is.na(res[,1]), "log2FC"] <- NA
      if(!is.numeric(predictor)){
        res$ordering <- NA
        res[!is.na(res[,coeff]) & res[,coeff] > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[1])
        res[!is.na(res[,coeff]) & res[,coeff] < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[coeff])
      }
    }
    
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- rownames(res)
    res$Method <- "Log Linear reg. 2 (llm2)"
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    return(res)
  }  
}
