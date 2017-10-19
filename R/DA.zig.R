#' MetageonomeSeq ZIG
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param by Column number or column name specifying which coefficient or contrast of the linear model is of interest (only for categorical predictors). Default 2
#' @param eff Filter features to have at least a "eff" quantile or number of effective samples.
#' @param allResults If TRUE will return raw results from the fitZig function
#' @param ... Additional arguments for the fitZig function
#' @export

DA.zig <- function(data, predictor, paired = NULL, covars = NULL, p.adj = "fdr", by = 2, eff = 0, allResults = FALSE, ...){
  
  suppressMessages(library(metagenomeSeq))
  
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

  count_table <- as.data.frame.matrix(count_table)
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  if(!is.null(paired)){
    if(is.null(covars)){
      mod <- model.matrix(~ predictor+paired)
    } else {
      if(class(data) == "phyloseq"){
        mod <- model.matrix(as.formula(paste("~ predictor+paired+",paste(covars, collapse="+"),sep = "")))
      } else {
        mod <- model.matrix(as.formula(paste("~ predictor+paired+",paste(names(covars), collapse="+"),sep = "")))
      }
    }
  } else {
    if(is.null(covars)){
      mod <- model.matrix(~ predictor)
    } else {
      if(class(data) == "phyloseq"){
        mod <- model.matrix(as.formula(paste("~ predictor+",paste(covars, collapse="+"),sep = "")))
      } else {
        mod <- model.matrix(as.formula(paste("~ predictor+",paste(names(covars), collapse="+"),sep = "")))
      }
    }
  }
  
  mgsfit <- fitZig(obj=mgsdata,mod=mod)
  if(is.numeric(predictor[1])){
    temp_table <- MRtable(mgsfit, number=nrow(count_table), by = by, coef = 1:2, eff = eff)
  } else {
    temp_table <- MRtable(mgsfit, number=nrow(count_table), by = by, coef = c(1:length(levels(as.factor(predictor)))), eff = eff)
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
  
  if(allResults) return(fitZig) else return(temp_table)
}


