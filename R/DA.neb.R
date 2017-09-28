#' Negative binomial glm
#'
#' With log(librarySize) as offset.
#' Mixed-effect model is used when a paired argument is included, with the paired variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param out.anova If TRUE will output results and p-values from anova (drop1 if paired != NULL). If false will output results for 2. level of the predictor.
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the glm.nb/glmer.nb function
#' @param ... Additional arguments for the glm.nb/glmer.nb functions
#' @import MASS
#' @importFrom lme4 glmer.nb glmer
#' @export

DA.neb <- function(data, predictor, paired = NULL, covars = NULL, relative = TRUE, out.anova = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
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
  
  if(relative) libSize <- colSums(count_table) else libSize <- 1
  count_table <- as.data.frame.matrix(count_table)
  
  if(is.null(paired)){
    if(is.null(covars)){
      negbin <- function(x){
        fit <- NULL
        tryCatch(
          fit <- MASS::glm.nb(x ~ predictor + offset(log(libSize)),...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(coef(summary(fit))) > 1) {
            coef(summary(fit))[2,]
          } else NA
        } else NA 
      }
    } else {
      if(class(data) == "phyloseq"){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),...), 
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
            fit <- MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(coef(summary(fit))) > 1) {
              coef(summary(fit))[2,]
            } else NA
          } else NA 
        }
      }
    }
  } else {
    if(is.null(covars)){
      negbin <- function(x){
        fit <- NULL
        tryCatch(
          fit <- lme4::glmer.nb(x ~ predictor + offset(log(libSize)) + (1|paired), ...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(coef(summary(fit))) > 1) {
            coef(summary(fit))[2,]
          } else NA
        } else NA 
      } 
    } else {
      if(class(data) == "phyloseq"){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(covars, collapse="+"),sep = "")), ...), 
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
            fit <- lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(names(covars), collapse="+"),sep = "")), ...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(coef(summary(fit))) > 1) {
              coef(summary(fit))[2,]
            } else NA
          } else NA 
        } 
      }
    }
  }
  
  if(out.anova){
    if(is.null(paired)){
      if(is.null(covars)){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- anova(MASS::glm.nb(x ~ predictor + offset(log(libSize)),...),test="Chisq")[2,], 
            error = function(x) fit <- NULL)
        }
      } else {
        if(class(data) == "phyloseq"){
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- anova(MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),...),test="Chisq")[2,], 
              error = function(x) fit <- NULL)
          }
        } else {
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- anova(MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),...),test="Chisq")[2,], 
              error = function(x) fit <- NULL)
          }
        }
      }
    } else {
      if(is.null(covars)){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- drop1(lme4::glmer.nb(x ~ predictor + offset(log(libSize)) + (1|paired), ...), test = "Chisq")[2,], 
            error = function(x) fit <- NULL)
        } 
      } else {
        if(class(data) == "phyloseq"){
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- drop1(lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(covars, collapse="+"),sep = "")), ...), test = "Chisq")[2,], 
              error = function(x) fit <- NULL)
          } 
        } else {
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- drop1(lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(names(covars), collapse="+"),sep = "")), ...), test = "Chisq")[2,], 
              error = function(x) fit <- NULL)
          } 
        }
      }
    }
  }
  
  if(out.anova){
    if(is.null(paired)){
      res <- as.data.frame(do.call(rbind,apply(count_table,1,negbin)))
      colnames(res) <- c("Df","Deviance","Resid. Df","Resid. Dev","pval")
    } else {
      res <- as.data.frame(do.call(rbind,apply(count_table,1,negbin)))
      colnames(res) <- c("Df","AIC","LRT","pval")
    }
  } else {
    res <- as.data.frame(t(as.data.frame(apply(count_table,1,negbin))))
    colnames(res) <- c("Estimate","Std.Error","t-value","pval")
  }
  
  
  if(nrow(res) == 1){
    res <- data.frame(Estimate = rep(NA,nrow(count_table)), Std.Error = rep(NA,nrow(count_table)), z.value = rep(NA,nrow(count_table)), pval = rep(NA,nrow(count_table)))
    rownames(res) <- rownames(count_table)                                                                                                           
  } 
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Negbinom GLM (neb)"
  
  if(nrow(res) > 1){
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  }

  if(allResults){
    if(is.null(paired)){
      if(is.null(covars)){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- MASS::glm.nb(x ~ predictor + offset(log(libSize)),...), 
            error = function(x) fit <- NULL)
        }
      } else {
        if(class(data) == "phyloseq"){
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),...), 
              error = function(x) fit <- NULL)
          }
        } else {
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- MASS::glm.nb(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),...), 
              error = function(x) fit <- NULL)
          }
        }
      }
    } else {
      if(is.null(covars)){
        negbin <- function(x){
          fit <- NULL
          tryCatch(
            fit <- lme4::glmer.nb(x ~ predictor + offset(log(libSize)) + (1|paired), ...), 
            error = function(x) fit <- NULL)
        } 
      } else {
        if(class(data) == "phyloseq"){
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(covars, collapse="+"),sep = "")), ...), 
              error = function(x) fit <- NULL)
          } 
        } else {
          negbin <- function(x){
            fit <- NULL
            tryCatch(
              fit <- lme4::glmer.nb(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(names(covars), collapse="+"),sep = "")), ...), 
              error = function(x) fit <- NULL)
          } 
        }
      }
    }
    return(apply(count_table,1,negbin))
  } else {
    return(res)
  }
  
}
