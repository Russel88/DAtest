#' EdgeR quasi-likelihood - TMM normalization
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the glmQLFTest function
#' @param ... Additional arguments for the calcNormFactors, estimateDisp, glmQLFit and glmQLFTest functions
#' @export

DA.erq <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", allResults = FALSE, ...){
  
  suppressMessages(library(edgeR))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- unlist(sample_data(data)[,predictor])
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    if(!is.null(covars)){
      for(i in 1:length(covars)){
        assign(covars[i], unlist(sample_data(data)[,covars[i]]))
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
  
  DA.erq.args <- list(...)
  calcNormFactors.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(calcNormFactors))]
  estimateDisp.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(estimateDisp))]
  glmQLFit.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFit))]
  glmQLFTest.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFTest))]
  
  
  count_table <- as.data.frame(count_table)
  y <- DGEList(counts=count_table,genes = data.frame(Feature = row.names(count_table)))
  y <- do.call(edgeR::calcNormFactors, c(list(y, method = "TMM"),calcNormFactors.args))
  if(!is.null(paired)){
    if(is.null(covars)){
      design <- model.matrix(~ predictor+paired)
    } else {
      if(class(data) == "phyloseq"){
        design <- model.matrix(as.formula(paste("~ predictor+paired+",paste(covars, collapse="+"),sep = "")))
      } else {
        design <- model.matrix(as.formula(paste("~ predictor+paired+",paste(names(covars), collapse="+"),sep = "")))
      }
    }
  } else {
    if(is.null(covars)){
      design <- model.matrix(~ predictor)
    } else {
      if(class(data) == "phyloseq"){
        design <- model.matrix(as.formula(paste("~ predictor+",paste(covars, collapse="+"),sep = "")))
      } else {
        design <- model.matrix(as.formula(paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")))
      }
    }
  }
  y <- do.call(estimateDisp,c(list(y,design),estimateDisp.args))
  fit <- do.call(glmQLFit,c(list(y,design),glmQLFit.args))
  
  if(is.numeric(predictor[1])){
    qlf <- do.call(glmQLFTest,c(list(fit,coef=2),glmQLFTest.args))
    ta <- qlf$table
    colnames(ta)[4] <- "pval"
  } else {
    qlf <- do.call(glmQLFTest,c(list(fit,coef=2:length(levels(as.factor(predictor)))),glmQLFTest.args))
    ta <- qlf$table
    colnames(ta)[(2+length(levels(as.factor(predictor))))] <- "pval"
  }

  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR qll - TMM (erq)"
  
  if(class(data) == "phyloseq") ta <- add.tax.DA(data, ta)
  
  if(allResults) return(qlf) else return(ta)
  
}


