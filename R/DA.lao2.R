#' ANOVA
#' 
#' With log transformation of relative abundances.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param covars Either a named list with covariables, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param delta Numeric. Pseudocount for the log transformation. Default 0.001
#' @param ... Additional arguments for the aov functions
#' @export

DA.lao2 <- function(data, predictor, covars = NULL, p.adj = "fdr", delta = 0.001, ...){
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
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
  
  if(is.null(covars)){
    form <- paste("x ~ predictor")
  } else {
    if(class(data) == "phyloseq"){
      form <- paste("x ~ ",paste(covars, collapse="+"),"+ predictor",sep = "")
    } else {
      form <- paste("x ~ ",paste(names(covars), collapse="+"),"+ predictor",sep = "")
    }
  }
  
  ao <- function(x){
    tryCatch(as.numeric(summary(aov(as.formula(form), ...))[[1]][(length(covars)+1),5]), error = function(e){NA}) 
  }
  
  count.rel <- apply(count_table,2,function(x) x/sum(x))
  count.rel <- log(count.rel + delta)
  
  res <- data.frame(pval = apply(count.rel,1,ao))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
    res$Feature <- rownames(res)
  res$Method <- "Log ANOVA 2"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  return(res)
}


