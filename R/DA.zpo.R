#' Zero inflated Poisson glm
#'
#' With log(librarySize) as offset.
#' Mixed-effect model is used when a paired argument is included, with the paired variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param out.anova If TRUE will output results and p-values from drop1. If false will output results for 2. level of the predictor.
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the zeroinfl function
#' @param ... Additional arguments for the zeroinfl function
#' @import pscl
#' @export

DA.zpo <- function(data, predictor, covars = NULL, out.anova = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
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
  
  libSize <- colSums(count_table)
  count_table <- as.data.frame.matrix(count_table)
  
  if(is.null(covars)){
      pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- pscl::zeroinfl(x ~ predictor + offset(log(libSize)),dist="poisson",...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(summary(fit)$coefficients$count) > 1) {
            summary(fit)$coefficients$count[2,]
          } else NA
        } else NA 
      }
    } else {
      if(class(data) == "phyloseq"){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),dist="poisson",...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(summary(fit)$coefficients$count) > 1) {
              summary(fit)$coefficients$count[2,]
            } else NA
          } else NA 
        }
      } else {
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),dist="poisson",...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(summary(fit)$coefficients$count) > 1) {
              summary(fit)$coefficients$count[2,]
            } else NA
          } else NA 
        }
      }
    }

  if(out.anova){
    if(is.null(covars)){
      pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- drop1(pscl::zeroinfl(x ~ predictor + offset(log(libSize)),dist="poisson",...),test="Chisq")[2,], 
          error = function(x) fit <- NULL)
      }
    } else {
      if(class(data) == "phyloseq"){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- drop1(pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),dist="poisson",...),test="Chisq")[2,], 
            error = function(x) fit <- NULL)
        }
      } else {
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- drop1(pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),dist="poisson",...),test="Chisq")[2,], 
            error = function(x) fit <- NULL)
        }
      }
    }
    res <- as.data.frame(do.call(rbind,apply(count_table,1,pois)))
    colnames(res) <- c("Df","AIC","LRT","pval")
  } else {
    res <- as.data.frame(t(as.data.frame(apply(count_table,1,pois))))
    colnames(res) <- c("Estimate","Std.Error","z value","pval")
  }
  
  
  if(nrow(res) == 1){
    res <- data.frame(Estimate = rep(NA,nrow(count_table)), Std.Error = rep(NA,nrow(count_table)), z.value = rep(NA,nrow(count_table)), pval = rep(NA,nrow(count_table)))
    rownames(res) <- rownames(count_table) 
    colnames(res) <- c("Estimate","Std.Error","z value","pval")
  } 
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "ZI-Poisson GLM (zpo)"
  
  if(nrow(res) > 1){
    if(class(data) == "phyloseq"){
      if(!is.null(tax_table(data, errorIfNULL = FALSE))){
        tax <- tax_table(data)
        res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
        rownames(res) <- NULL
      } 
    }
  }

  if(allResults){
    if(is.null(covars)){
      pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- pscl::zeroinfl(x ~ predictor + offset(log(libSize)),dist="poisson",...), 
          error = function(x) fit <- NULL)
      }
    } else {
      if(class(data) == "phyloseq"){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),dist="poisson",...), 
            error = function(x) fit <- NULL)
        }
      } else {
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- pscl::zeroinfl(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),dist="poisson",...), 
            error = function(x) fit <- NULL)
        }
      }
    }
    
    return(apply(count_table,1,pois))
  } else {
    return(res)
  }
  
}
