#' Poisson glm
#'
#' With librarySize as offset.
#' Mixed-effect model is used when a paired argument is included, with the paired variable as a random intercept.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the glm/glmer functions
#' @importFrom lme4 glmer
#' @export

DA.poi <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", ...){
 
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
  
  libSize <- colSums(count_table)
  count_table <- as.data.frame.matrix(count_table)
  
  if(is.null(paired)){
    if(is.null(covars)){
      pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- glm(x ~ predictor + offset(log(libSize)),family="poisson",...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(coef(summary(fit))) > 1) {
            coef(summary(fit))[2,]
          } else NA
        } else NA 
      }
    } else {
      if(class(data) == "phyloseq"){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- glm(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(covars, collapse="+"),sep = "")),family="poisson",...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(coef(summary(fit))) > 1) {
              coef(summary(fit))[2,]
            } else NA
          } else NA 
        }
      } else {
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- glm(as.formula(paste("x ~ predictor+offset(log(libSize))+",paste(names(covars), collapse="+"),sep = "")),family="poisson",...), 
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
      pois <- function(x){
        fit <- NULL
        tryCatch(
          fit <- lme4::glmer(x ~ predictor + offset(log(log(libSize))) + (1|paired),family="poisson", ...), 
          error = function(x) fit <- NULL)
        if(!is.null(fit)) {
          if(nrow(coef(summary(fit))) > 1) {
            coef(summary(fit))[2,]
          } else NA
        } else NA 
      } 
    } else {
      if(class(data) == "phyloseq"){
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- lme4::glmer(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(covars, collapse="+"),sep = "")),family="poisson", ...), 
            error = function(x) fit <- NULL)
          if(!is.null(fit)) {
            if(nrow(coef(summary(fit))) > 1) {
              coef(summary(fit))[2,]
            } else NA
          } else NA 
        } 
      } else {
        pois <- function(x){
          fit <- NULL
          tryCatch(
            fit <- lme4::glmer(as.formula(paste("x ~ predictor+offset(log(libSize)) + (1|paired)+",paste(names(covars), collapse="+"),sep = "")),family="poisson", ...), 
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
  
  res <- as.data.frame(t(as.data.frame(apply(count_table,1,pois))))
  if(nrow(res) == 1){
    res <- data.frame(Estimate = rep(NA,nrow(count_table)), Std.Error = rep(NA,nrow(count_table)), z.value = rep(NA,nrow(count_table)), pval = rep(NA,nrow(count_table)))
    rownames(res) <- rownames(count_table)                                                                                                           
  } 
  colnames(res) <- c("Estimate","Std.Error","z value","pval")
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- rownames(res)
  res$Method <- "Poisson GLM"
  
  if(nrow(res) > 1){
    if(class(data) == "phyloseq"){
      if(!is.null(tax_table(data, errorIfNULL = FALSE))){
        tax <- tax_table(data)
        res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
        rownames(res) <- NULL
      } 
    }
  }

  return(res)
  
}
