#' Mvabund
#'
#' Implementation of mvabund manyglm for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param relative Logical. Whether \code{log(librarySize)} should be used as offset. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details. Alternatively, "mva" uses mvabunds adjusted p-values
#' @param coeff Integer. The p-value and log2FoldChange will be associated with this coefficient. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param coeff.ref Integer. Reference level of the \code{predictor}. Will only affect the log2FC and ordering columns on the output. Default the intercept, = 1 
#' @param resamp Resample method for estimating p-values. Passed to \code{summary.manyglm}. Default "perm.resid"
#' @param allResults If TRUE will return raw results from the \code{mvabund} function
#' @param ... Additional arguments for the \code{manyglm} and \code{summary.manyglm} functions
#' @export

DA.mva <- function(data, predictor, paired = NULL, covars = NULL, relative = TRUE, p.adj = "fdr", coeff = 2, coeff.ref = 1, resamp = "perm.resid", allResults = FALSE, ...){
  
  suppressMessages(library(mvabund))
  
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

  # Extract arguments
  DA.mva.args <- list(...)
  manyglm.args <- DA.mva.args[names(DA.mva.args) %in% names(formals(manyglm))]
  summary.manyglm.args <- DA.mva.args[names(DA.mva.args) %in% names(formals(summary.manyglm))]
  
  mva.table <- mvabund(t(count_table))

  # Define model
  if(is.null(covars)){
    if(is.null(paired)){
      form <- as.formula(mva.table ~ predictor)
    } else {
      form <- as.formula(mva.table ~ predictor + paired)
    }
  } else {
    if(is.null(paired)){
      form <- as.formula(paste("mva.table ~ predictor+",paste(names(covars), collapse="+"),sep = ""))
    } else {
      form <- as.formula(paste("mva.table ~ predictor+paired+",paste(names(covars), collapse="+"),sep = ""))
    }
  }
  
  # Offset
  if(relative){
    libSize <- colSums(count_table)
    form <- as.formula(paste(as.list(form)[[2]],"~",as.list(form[3:length(as.list(form))]),"+offset(log(libSize))"))
  }
  
  # Fit model
  mod <- do.call(manyglm,c(list(form, ...),manyglm.args))
  
  if(allResults){
    return(mod)
  } else {
    pvals.unadj <- do.call(summary,c(list(mod, p.uni="unadjusted", resamp = resamp, ...),summary.manyglm.args))
    
    # Extract results
    res <- as.data.frame(cbind(t(mod$coefficients),pvals.unadj$uni.p[,2]))
    colnames(res)[ncol(res)] <- c("pval")
    if(p.adj == "mva"){
      pvals.adj <- do.call(summary,c(list(mod, p.uni="adjusted", resamp = resamp, ...),summary.manyglm.args))
      res$pval.adj <- pvals.adj$uni.p[,2]
    } else {
      res$pval.adj <- p.adjust(res$pval, method = p.adj)
    } 
    
    res$log2FC <- log2(exp(res[,coeff.ref]+res[,coeff]) / exp(res[,coeff.ref]))
    
    if(!is.numeric(predictor)){
      res$ordering <- NA
      res[!is.na(res[,coeff]) & res[,coeff] > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[coeff.ref])
      res[!is.na(res[,coeff]) & res[,coeff] < 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff.ref],">",levels(as.factor(predictor))[coeff])
    }
    
    res$Feature <- rownames(pvals.unadj$uni.test)
    res$Method <- "mvabund (mva)"
    
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    
    return(res)
  }
  
}


