#' LIMMA with \code{voom}
#' 
#' Limma is run after variance stabilization with \code{voom} using \code{edgeR} normalized library sizes
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE will output results from F-tests, if FALSE t-statistic results from 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The p-value and log2FoldChange will be associated with this coefficient. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{eBayes} function
#' @param ... Additional arguments for the \code{voom}, \code{eBayes} and \code{lmFit} functions
#' @export

DA.vli <- function(data, predictor, paired = NULL, covars = NULL, out.all = NULL, p.adj = "fdr", coeff = 2, allResults = FALSE, ...){
  
  suppressMessages(library(limma))
  if(!is.null(paired)) suppressMessages(library(statmod))

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
  limma.args <- list(...)
  lmFit.args <- limma.args[names(limma.args) %in% names(formals(lmFit))]
  eBayes.args <- limma.args[names(limma.args) %in% names(formals(eBayes))]
  voom.args <- limma.args[names(limma.args) %in% names(formals(voom))]
  
  count_table <- as.data.frame(count_table)
  
  # Define model
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
  
  # Variance stabilization
  nf <- edgeR::calcNormFactors(count_table)
  voo <- do.call(voom,c(list(count_table, design, lib.size=colSums(count_table)*nf),voom.args))
  
  # Fit linear model
  if(is.null(paired)){
    fit <- do.call(lmFit,c(list(voo, design),lmFit.args))
  } else {
    dupcor <-  duplicateCorrelation(count_table, design, block = paired)
    fit <- do.call(lmFit,c(list(voo, design, block = paired, correlation = dupcor$cor),lmFit.args))
  }
  
  # Fit empirical bayes
  fit.eb <- do.call(eBayes, c(list(fit),eBayes.args))
  
  # Extract results
  if(is.numeric(predictor[1])){
      res <- topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = 2)
      colnames(res)[4:5] <- c("pval","pval.adj")
  } else {
    if(out.all){
      res <- topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = 2:length(levels(as.factor(predictor))))
      colnames(res)[length(levels(as.factor(predictor)))+2:3] <- c("pval","pval.adj")
    } else {
      res <- topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = coeff)
      colnames(res)[4:5] <- c("pval","pval.adj")
      res$ordering <- NA
      res[!is.na(res$logFC) & res$logFC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[1])
      res[!is.na(res$logFC) & res$logFC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[coeff])
    }
  }
  
  res$Feature <- rownames(res)
  res$Method <- "LIMMA voom (vli)"
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)

  if(allResults) return(fit.eb) else return(res)  
}

