#' EdgeR quasi-likelihood
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the calcNormFactors, estimateDisp, glmQLFit and glmQLFTest functions
#' @export

DA.erq <- function(data, predictor, paired = NULL, p.adj = "fdr", ...){
  
  library(edgeR)
  
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
  
  DA.erq.args <- list(...)
  calcNormFactors.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(calcNormFactors))]
  estimateDisp.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(estimateDisp))]
  glmQLFit.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFit))]
  glmQLFTest.args <- DA.erq.args[names(DA.erq.args) %in% names(formals(glmQLFTest))]
  
  
  count_table <- as.data.frame(count_table)
  y <- DGEList(counts=count_table,genes = data.frame(Feature = row.names(count_table)))
  y <- do.call(calcNormFactors, c(list(y),calcNormFactors.args))
  if(is.null(paired)){
    design <- model.matrix(~predictor)
  } else {
    design <- model.matrix(~predictor + paired)
  }
  y <- do.call(estimateDisp,c(list(y,design),estimateDisp.args))
  fit <- do.call(glmQLFit,c(list(y,design),glmQLFit.args))
  qlf <- do.call(glmQLFTest,c(list(fit,coef=2),glmQLFTest.args))
  ta <- qlf$table
  colnames(ta)[4] <- "pval"
  ta$pval.adj <- p.adjust(ta$pval, method = p.adj)
  ta$Feature <- rownames(ta)
  ta$Method <- "EdgeR qll"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      ta <- merge(ta, tax, by.x = "Feature", by.y = "row.names")
      rownames(ta) <- NULL
    } 
  }
  
  return(ta)
  
}


