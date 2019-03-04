#' MetagenomeSeq ZIG
#'
#' Implementation of Metagenome zero-inflated gaussian model for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param by Column number or column name specifying which coefficient or contrast of the linear model is of interest (only for categorical \code{predictor}). Default 2
#' @param eff Filter features to have at least a \code{eff} quantile or number of effective samples, passed to \code{MRtable}
#' @param allResults If TRUE will return raw results from the \code{fitZig} function
#' @param ... Additional arguments for the \code{fitZig} function
#' @return A data.frame with with results.
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(1000, size = 0.1, mu = 500), nrow = 100, ncol = 10)
#' rownames(mat) <- 1:100
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running MetagenomeSeq Zero-inflated Gaussian
#' res <- DA.zig(data = mat, predictor = pred)
#' @export

DA.zig <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", by = 2, eff = 0.5, allResults = FALSE, ...){
  
  ok <- tryCatch({
    loadNamespace("metagenomeSeq")
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
    
    # Collect data and normalize
    count_table <- as.data.frame.matrix(count_table)
    mgsdata <- metagenomeSeq::newMRexperiment(counts = count_table)
    mgsp <- metagenomeSeq::cumNormStat(mgsdata)
    mgsdata <- metagenomeSeq::cumNorm(mgsdata, mgsp)
    
    # Define model
    if(is.null(covars)){
      mod <- model.matrix(~ predictor)
    } else {
      mod <- model.matrix(as.formula(paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")))
    }
    
    # Fit model
    if(is.null(paired)){
      mgsfit <- metagenomeSeq::fitZig(obj=mgsdata,mod=mod,...)
    } else {
      mgsfit <- metagenomeSeq::fitZig(obj=mgsdata,mod=mod,...,useMixedModel=TRUE,block=paired)
    }
    
    # Extract results
    if(is.numeric(predictor)){
      temp_table <- metagenomeSeq::MRtable(mgsfit, number=nrow(count_table), by = by, coef = 1:2, eff = eff)
      colnames(temp_table)[6] <- "logFC"
      temp_table <- temp_table[,-ncol(temp_table)]
    } else {
      temp_table <- metagenomeSeq::MRtable(mgsfit, number=nrow(count_table), by = by, coef = c(seq_along(levels(as.factor(predictor)))), eff = eff)
      temp_table <- temp_table[,-ncol(temp_table)]
      if(length(levels(as.factor(predictor))) == 2){
        colnames(temp_table)[6] <- "logFC"
        temp_table$ordering <- NA
        temp_table[!is.na(temp_table$logFC) & temp_table$logFC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[by],">",levels(as.factor(predictor))[1])
        temp_table[!is.na(temp_table$logFC) & temp_table$logFC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[by])
      }
    }
    temp_table <- temp_table[!is.na(row.names(temp_table)),]
    
    # Pvalue have different naming depending on package version
    if("pvalues" %in% names(temp_table)){
      colnames(temp_table)[which(names(temp_table) == "pvalues")] <- "pval"
    } 
    if("pValue" %in% names(temp_table)){
      colnames(temp_table)[which(names(temp_table) == "pValue")] <- "pval"
    }
    temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
    temp_table$Feature <- rownames(temp_table)
    temp_table$Method <- "MgSeq ZIG (zig)"
    
    if(class(data) == "phyloseq") temp_table <- add.tax.DA(data, temp_table)
    
    if(allResults) return(mgsfit) else return(temp_table)
    
  } else {
    stop("metagenomeSeq package required")
  }
  
}


