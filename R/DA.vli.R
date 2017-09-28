#' LIMMA voom
#'
#' Some implementation is borrowed from:
#' http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param out.anova If TRUE will output results from F-tests, if FALSE t-statistic results from 2. level of the predictor.
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the eBayes function
#' @param ... Additional arguments for the voom, eBayes and lmFit functions
#' @import statmod
#' @export

DA.vli <- function(data, predictor, paired = NULL, covars = NULL, out.anova = TRUE, p.adj = "fdr", allResults = FALSE, ...){
  
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
  
  limma.args <- list(...)
  lmFit.args <- limma.args[names(limma.args) %in% names(formals(lmFit))]
  eBayes.args <- limma.args[names(limma.args) %in% names(formals(eBayes))]
  voom.args <- limma.args[names(limma.args) %in% names(formals(voom))]
  
  if(is.null(covars)){
    form <- paste("~ predictor")
  } else {
    if(class(data) == "phyloseq"){
      form <- paste("~ predictor+",paste(covars, collapse="+"),sep = "")
    } else {
      form <- paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")
    }
  }
  
  design <- model.matrix(as.formula(form))
  nf <- edgeR::calcNormFactors(count_table)
  voo <- do.call(voom,c(list(count_table, design, lib.size=colSums(count_table)*nf),voom.args))
  n <- dim(count_table)[1]
  if(is.null(paired)){
    fit <- do.call(lmFit,c(list(voo, design),lmFit.args))
  } else {
    dupcor <-  duplicateCorrelation(count_table, design, block = paired)
    fit <- do.call(lmFit,c(list(voo, design, block = paired, correlation = dupcor$cor),lmFit.args))
  }
  fit.eb <- do.call(eBayes, c(list(fit),eBayes.args))
  
  if(out.anova){
    if(is.numeric(predictor[1])){
      res <- topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = 2)
      colnames(res)[4:5] <- c("pval","pval.adj")
    } else {
      res <- topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = 2:length(levels(as.factor(predictor))))
      colnames(res)[length(levels(as.factor(predictor)))+2:3] <- c("pval","pval.adj")
    }
  } else {
    Estimate <- fit.eb$coefficients
    df.residual <- fit.eb$df.residual
    df.prior <- rep(fit.eb$df.prior, n)
    s2.prior <- rep(fit.eb$s2.prior, n)
    s2 <- (fit.eb$sigma)^2
    s2.post <- fit.eb$s2.post
    stat <- fit.eb$t[,2]
    pval <- fit.eb$p.value[,2]
    pval.adj <- p.adjust(pval, method = p.adj)
    res <- data.frame(Estimate, stat, pval, pval.adj, df.residual, df.prior, s2.prior, s2, s2.post)
  }
  
  res$Feature <- rownames(res)
  res$Method <- "LIMMA voom (vli)"
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)

  if(allResults) return(fit.eb) else return(res)  
}

