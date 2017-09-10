#' LIMMA
#' 
#' With log transformation of relative abundances.
#' Some implementation is borrowed from:
#' http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for log transformation. Default 0.001
#' @param ... Additional arguments for the eBayes and lmFit functions
#' @export

DA.lli2 <- function(data, predictor, paired = NULL, p.adj = "fdr", delta = 0.001, correlation = 0.75, ...){
  
  library(limma)
  library(statmod)
  
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
  
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  
  count.rel <- log(count.rel + delta)
    
  limma.args <- list(...)
  lmFit.args <- limma.args[names(limma.args) %in% names(formals(lmFit))]
  eBayes.args <- limma.args[names(limma.args) %in% names(formals(eBayes))]
  
  design <- model.matrix(~predictor)
  n <- dim(count.rel)[1]
  if(is.null(paired)){
    fit <- do.call(lmFit,c(list(count.rel, design),lmFit.args))
  } else {
    dupcor <-  duplicateCorrelation(count.rel, design, block = paired)
    fit <- do.call(lmFit,c(list(count.rel, design, block = paired, correlation = dupcor$cor),lmFit.args))
  }
  fit.eb <- do.call(eBayes, c(list(fit),eBayes.args))
  Estimate <- fit.eb$coefficients
  df.residual <- fit.eb$df.residual
  df.prior <- rep(fit.eb$df.prior, n)
  s2.prior <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.stat <- fit.eb$t[, 2]
  pval <- fit.eb$p.value[, 2]
  pval.adj <- p.adjust(pval, method = p.adj)
  res <- data.frame(Estimate, t.stat, pval, pval.adj, df.residual, df.prior, s2.prior, s2, s2.post)
  res$Feature <- rownames(res)
  res$Method <- "Log LIMMA 2"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }

  return(res)  
}

