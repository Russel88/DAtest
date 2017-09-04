#' EdgeR quasi-likelihood
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param outcome The outcome of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the glmQLFit function
#' @export

DA.erq <- function(data, outcome, paired = NULL, p.adj = "fdr", ...){
  
  library(edgeR)
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(outcome) > 1 | length(paired) > 1) stop("When data is a phyloseq object outcome and paired should only contain the name of the variables in sample_data")
    if(!outcome %in% sample_variables(data)) stop(paste(outcome,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    outcome <- suppressWarnings(as.matrix(sample_data(data)[,outcome]))
    if(!is.null(paired)) paired <- suppressWarnings(as.matrix(sample_data(data)[,paired]))
  } else {
    count_table <- data
  }
  
  otu_table <- as.data.frame(count_table)
  y <- DGEList(counts=count_table,genes = data.frame(Feature = row.names(count_table)))
  y <- calcNormFactors(y)
  if(is.null(paired)){
    design <- model.matrix(~outcome)
  } else {
    design <- model.matrix(~outcome + paired)
  }
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design, ...)
  if(is.numeric(outcome)){
    qlf <- glmQLFTest(fit,coef=2)
    ta <- qlf$table
    colnames(ta)[4] <- "pval"
  } else {
    qlf <- glmQLFTest(fit,coef=seq(2,length(levels(as.factor(outcome)))))
    ta <- qlf$table
    colnames(ta)[(length(levels(as.factor(outcome)))+2)] <- "pval"
  }

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


