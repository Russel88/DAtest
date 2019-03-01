#' DESeq2 with manual geometric means
#'
#' Implementation of DESeq2 for \code{DAtest}
#' Manual geometric means calculated to avoid errors, see https://github.com/joey711/phyloseq/issues/387
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE, will run "LRT" which will produce one p-value for the \code{predictor}. If FALSE will run "Wald" test and will output p-value from one level of the predictor specified by \code{coeff}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The log2FoldChange (and p-value if test="Wald") will be associated with this coefficient. This coefficient is by default compared to the intercept (1. level of \code{predictor}), change this with \code{coeff.ref}. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param coeff.ref Integer. Reference level of the \code{predictor}. Default the intercept, = 1 
#' @param allResults If TRUE will return raw results from the \code{DESeq} function
#' @param ... Additional arguments for the \code{DESeq} function
#' @return A data.frame with with results.
#' @export

DA.ds2 <- function(data, predictor, paired = NULL, covars = NULL, out.all = NULL, p.adj = "fdr", coeff = 2, coeff.ref = 1, allResults = FALSE, ...){
  
  ok <- tryCatch({
    loadNamespace("DESeq2")
    TRUE
  }, error=function(...) FALSE)
  
  if (ok){
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
    predictor <- as.factor(predictor)
    
    if(coeff == coeff.ref) stop("coeff and coeff.ref cannot be the same")
    if(!coeff %in% seq_along(unique(predictor)) | !coeff.ref %in% seq_along(unique(predictor))) stop(paste("coeff and coeff.ref should be integers between 1 and",length(unique(predictor))))
    
    # out.all
    if(is.null(out.all)){
      if(length(unique(predictor)) == 2) out.all <- FALSE
      if(length(unique(predictor)) > 2) out.all <- TRUE
      if(is.numeric(predictor)) out.all <- FALSE
    }
    
    # Collect data
    if(is.null(paired)){
      if(is.null(covars)){
        predictordf <- data.frame(predictor = factor(predictor))
        row.names(predictordf) <- colnames(count_table)
        x <- DESeq2::DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = predictordf , design = ~ predictor)
      } else {
        predictordf <- as.data.frame(c(list(predictor = factor(predictor)),covars))
        row.names(predictordf) <- colnames(count_table)
        x <- DESeq2::DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = predictordf , design = as.formula(paste("~",paste(names(covars), collapse="+"),"+predictor",sep = "")))
      }
    } else {
      if(is.null(covars)){
        predictordf <- data.frame(predictor = factor(predictor),
                                  paired = factor(paired))
        row.names(predictordf) <- colnames(count_table)
        x <- DESeq2::DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = predictordf , design = ~ paired + predictor)
      } else {
        predictordf <- as.data.frame(c(list(predictor = factor(predictor),paired = factor(paired)),covars))
        row.names(predictordf) <- colnames(count_table)
        x <- DESeq2::DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = predictordf , design = as.formula(paste("~ paired+",paste(names(covars), collapse="+"),"+predictor",sep = "")))
      }
    }
    
    # Geometric means
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(DESeq2::counts(x), 1, gm_mean)
    x = DESeq2::estimateSizeFactors(x, geoMeans = geoMeans)
    
    # Run test
    if(out.all){
      if(is.null(paired)){
        if(is.null(covars)){
          x <- DESeq2::DESeq(x,test="LRT",reduced = ~1, ...)
        } else {
          x <- DESeq2::DESeq(x,test="LRT",reduced = as.formula(paste("~ ",paste(names(covars), collapse="+"),sep = "")), ...)
        } 
      } else {
        if(is.null(covars)){
          x <- DESeq2::DESeq(x,test="LRT",reduced = ~ paired, ...)
        } else {
          x <- DESeq2::DESeq(x,test="LRT",reduced = as.formula(paste("~ paired +",paste(names(covars), collapse="+"),sep = "")), ...)
        }
      }
      res <- as.data.frame(DESeq2::results(x, contrast = c("predictor",levels(predictor)[coeff],levels(predictor)[coeff.ref]))@listData)
      res$ordering <- NA
      res[!is.na(res$log2FoldChange) & res$log2FoldChange > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[coeff.ref])
      res[!is.na(res$log2FoldChange) & res$log2FoldChange < 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff.ref],">",levels(as.factor(predictor))[coeff])
    }
    if(!out.all){
      x <- DESeq2::DESeq(x,test="Wald", ...)
      if(paste0("predictor",levels(predictor)[coeff]) %in% DESeq2::resultsNames(x)){
        res <- as.data.frame(DESeq2::results(x, name = paste0("predictor",levels(predictor)[coeff]))@listData)
        res$ordering <- NA
        res[!is.na(res$log2FoldChange) & res$log2FoldChange > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">mean")
        res[!is.na(res$log2FoldChange) & res$log2FoldChange < 0,"ordering"] <- paste0("mean>",levels(as.factor(predictor))[coeff])
      } else {
        res <- as.data.frame(DESeq2::results(x, contrast = c("predictor",levels(predictor)[coeff],levels(predictor)[coeff.ref]))@listData)
        res$ordering <- NA
        res[!is.na(res$log2FoldChange) & res$log2FoldChange > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[coeff.ref])
        res[!is.na(res$log2FoldChange) & res$log2FoldChange < 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff.ref],">",levels(as.factor(predictor))[coeff])
      }
    }
    
    colnames(res)[5] <- "pval"
    res <- res[,-(ncol(res)-1)]
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- DESeq2::results(x)@rownames
    res$Method <- "DESeq2 man. geoMeans (ds2)"
    
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    
    if(allResults) return(x) else return(res)
  } else {
    stop("DESeq2 package required")
  }
    
}

