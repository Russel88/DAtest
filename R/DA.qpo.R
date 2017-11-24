#' Quasi-poisson glm
#'
#' Apply quasi-poisson generalized linear model for multiple features with one \code{predictor}
#' With \code{log(librarySize)} as offset if \code{relative=TRUE}.
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param relative Logical. Whether \code{log(librarySize)} should be used as offset. Default TRUE
#' @param out.all If TRUE will output results and p-values from \code{anova}. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The p-value and log2FoldChange will be associated with this coefficient. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{glm} function
#' @param ... Additional arguments for the \code{glm} functions
#' @export

DA.qpo <- function(data, predictor, covars = NULL, relative = TRUE, out.all = NULL, p.adj = "fdr", coeff = 2, allResults = FALSE, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired = NULL, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
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
  
  # Library sizes
  if(relative) libSize <- colSums(count_table) else libSize <- rep(1,ncol(count_table))
  count_table <- as.data.frame.matrix(count_table)

  # Define functions
  if(is.null(covars)){
    pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- glm(x ~ predictor + offset(log(libSize)),family="quasipoisson",...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(coef(summary(fit))) > 1) {
            pval <- coef(summary(fit))[coeff,4]
            ests <- coef(summary(fit))[,1]
            c(ests,pval)
          } else NA
        } else NA 
      }
  } else {
    pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- glm(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),family="quasipoisson",...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(coef(summary(fit))) > 1) {
              pval <- coef(summary(fit))[coeff,4]
              ests <- coef(summary(fit))[,1]
              c(ests,pval)
            } else NA
          } else NA 
        }
  }
  
  ## for out.all TRUE
  if(out.all){
    if(is.null(covars)){
      pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- glm(x ~ predictor + offset(log(libSize)),family = "quasipoisson",...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)){
            ests <- coef(summary(fit))[,1]
            ano <- anova(fit, test = "Chisq")[2,]
            c(ano,ests)
          }
        }
    } else {
      pois <- function(x){
            fit <- NULL
            tryCatch(
              fit <- glm(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),family = "quasipoisson",...), 
              error = function(x) fit <- NULL)
            if(!is.null(fit)){
              ests <- coef(summary(fit))[,1]
              ano <- anova(fit, test = "Chisq")[2,]
              c(ano,ests)
            }
          }
    }
  }

  # Run tests
  if(allResults){
      if(is.null(covars)){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- glm(x ~ predictor + offset(log(libSize)),family="quasipoisson",...), 
            error = function(x) fit <- NULL)
        }
      } else {
        pois <- function(x){
            fit <- NULL
            tryCatch(
              fit <- glm(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),family="quasipoisson",...), 
              error = function(x) fit <- NULL)
          }
      }
    return(apply(count_table,1,pois))
  } else {
    if(out.all){
      res <- as.data.frame(do.call(rbind,apply(count_table,1,pois)))
      colnames(res)[1:5] <- c("Df","Deviance","Resid. Df","Resid. Dev","pval")
      res <- as.data.frame(lapply(res, unlist))
    } else {
      res <- as.data.frame(t(as.data.frame(apply(count_table,1,pois))))
      colnames(res)[ncol(res)] <- "pval"
      res$log2FC <- log2(exp(res[,1]+res[,coeff]) / exp(res[,1]))
      if(!is.numeric(predictor)){
        res$ordering <- NA
        res[!is.na(res[,coeff]) & res[,coeff] > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[1])
        res[!is.na(res[,coeff]) & res[,coeff] < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[coeff])
      }
    }
    
    if(nrow(res) == 1){
      res <- data.frame(pval = rep(NA,nrow(count_table)))
      rownames(res) <- rownames(count_table)                                                                                                           
    } 
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- rownames(res)
    res$Method <- "Quasi-Poisson GLM (qpo)"
    
    if(nrow(res) > 1){
      if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    }
    return(res)
  }
  
}
