#' LIMMA
#' 
#' Limma with log transformation of relative abundances.
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE will output results from F-tests, if FALSE t-statistic results from 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param delta Numeric. Pseudocount for log transformation. Default 0.001
#' @param coeff Integer. The p-value and log2FoldChange will be associated with this coefficient. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{eBayes} function
#' @param ... Additional arguments for the \code{eBayes} and \code{lmFit} functions
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 100, ncol = 10)
#' rownames(mat) <- 1:100
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running limma
#' res <- DA.lli2(data = mat, predictor = pred)
#' @export

DA.lli2 <- function(data, predictor, paired = NULL, covars = NULL, out.all = NULL, p.adj = "fdr", delta = 0.001, coeff = 2, allResults = FALSE, ...){
  
  ok <- tryCatch({
    if(is.null(paired)){loadNamespace("limma")} else {loadNamespace(c("limma","statmod"))}
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
    if(!is.null(covars)){
      for(i in seq_along(covars)){
        assign(names(covars)[i], covars[[i]])
      }
    }
    
    # out.all
    if(is.null(out.all)){
      if(length(unique(predictor)) == 2) out.all <- FALSE
      if(length(unique(predictor)) > 2) out.all <- TRUE
      if(is.numeric(predictor)) out.all <- FALSE
    }
    
    # Relative abundance and log
    count.rel <- apply(count_table,2,function(x) x/sum(x))
    count.rel <- log(count.rel + delta)
    count.rel <- as.data.frame(count.rel)  
    
    # Arguments
    limma.args <- list(...)
    lmFit.args <- limma.args[names(limma.args) %in% names(formals(limma::lmFit))]
    eBayes.args <- limma.args[names(limma.args) %in% names(formals(limma::eBayes))]
    
    # Define design
    if(is.null(covars)){
      form <- paste("~ predictor")
    } else {
      form <- paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")
    }
    design <- model.matrix(as.formula(form))
    
    # Linear fit
    if(is.null(paired)){
      fit <- do.call(limma::lmFit,c(list(count.rel, design),lmFit.args))
    } else {
      dupcor <-  limma::duplicateCorrelation(count.rel, design, block = paired)
      fit <- do.call(limma::lmFit,c(list(count.rel, design, block = paired, correlation = dupcor$cor),lmFit.args))
    }
    
    # Empirical bayes
    fit.eb <- do.call(limma::eBayes, c(list(fit),eBayes.args))
    
    # Extract results
    if(is.numeric(predictor[1])){
      res <- limma::topTable(fit.eb, number = nrow(count.rel), adjust.method = p.adj, coef = 2)
      colnames(res)[4:5] <- c("pval","pval.adj")
    } else {
      if(out.all){
        res <- limma::topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = 2:length(levels(as.factor(predictor))))
        colnames(res)[length(levels(as.factor(predictor)))+2:3] <- c("pval","pval.adj")
      } else {
        res <- limma::topTable(fit.eb, number = nrow(count_table), adjust.method = p.adj, coef = coeff)
        colnames(res)[4:5] <- c("pval","pval.adj")
        res$ordering <- NA
        res[!is.na(res$logFC) & res$logFC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[coeff],">",levels(as.factor(predictor))[1])
        res[!is.na(res$logFC) & res$logFC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[coeff])
      }
    }
    
    res$Feature <- rownames(res)
    res$Method <- "Log LIMMA 2 (lli2)"
    
    if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
    
    if(allResults) return(fit.eb) else return(res) 
    
  } else {
    if(!is.null(paired)){
      stop("limma and statmod packages required")
    } else {
      stop("limma package required")
    }
  }
  
  
}

