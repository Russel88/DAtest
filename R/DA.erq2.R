#' EdgeR quasi-likelihood - RLE normalization
#' 
#' Implementation of edgeR quasi-likelihood test for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE will output one p-value for all levels of the predictor. If FALSE outputs p-value and logFC from one level of the \code{predictor} specified by \code{coeff}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The p-value and logFC will be associated with this coefficient when \code{out.all = FALSE}. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{glmQLFTest} function
#' @param ... Additional arguments for the \code{calcNormFactors}, \code{estimateDisp}, \code{glmQLFit} and \code{glmQLFTest} functions
#' @export

DA.erq2 <- function(data, predictor, paired = NULL, covars = NULL, out.all = NULL, p.adj = "fdr", coeff = 2, allResults = FALSE, ...){
  
  suppressMessages(library(edgeR))
  
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
  
  # Arguments
  DA.erq.args <- list(...)
  calcNormFactors.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(calcNormFactors))]
  estimateDisp.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(estimateDisp))]
  glmQLFit.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFit))]
  glmQLFTest.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFTest))]
  
  count_table <- as.data.frame(count_table)
  y <- DGEList(counts=count_table,genes = data.frame(Feature = row.names(count_table)))
  
  # Normalization
  y <- do.call(edgeR::calcNormFactors, c(list(y, method = "RLE"),calcNormFactors.args))
  
  # Define model
  if(!is.null(paired)){
    if(is.null(covars)){
      design <- model.matrix(~ predictor+paired)
    } else {
      design <- model.matrix(as.formula(paste("~ predictor+paired+",paste(names(covars), collapse="+"),sep = "")))
    }
  } else {
    if(is.null(covars)){
      design <- model.matrix(~ predictor)
    } else {
      design <- model.matrix(as.formula(paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")))
    }
  }
  
  # Dispersion and fit 
  y <- do.call(estimateDisp,c(list(y,design),estimateDisp.args))
  fit <- do.call(glmQLFit,c(list(y,design),glmQLFit.args))

  # Run tests
  if(is.numeric(predictor[1])){
    qlf <- do.call(glmQLFTest,c(list(fit,coef=2),glmQLFTest.args))
    ta <- qlf$table
    colnames(ta)[4] <- "pval"
  } else {
    if(out.all){
      qlf <- do.call(glmQLFTest,c(list(fit,coef=2:length(levels(as.factor(predictor)))),glmQLFTest.args))
      ta <- qlf$table
      colnames(ta)[(2+length(levels(as.factor(predictor))))] <- "pval"
    } else {
      qlf <- do.call(glmQLFTest,c(list(fit,coef=coeff),glmQLFTest.args))
      ta <- qlf$table
      colnames(ta)[4] <- "pval"
      ta$ordering <- NA
      ta[!is.na(ta$logFC) & ta$logFC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[1])
      ta[!is.na(ta$logFC) & ta$logFC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[coeff])
    }
  }
  
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR qll - RLE (erq2)"
  
  if(class(data) == "phyloseq") ta <- add.tax.DA(data, ta)
  
  if(allResults) return(qlf) else return(ta)
  
}


