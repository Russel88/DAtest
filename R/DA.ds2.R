#' DESeq2
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8.
#' Manual geometric means calculated to avoid errors, see https://github.com/joey711/phyloseq/issues/387
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param outcome The outcome of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the DESeq function
#' @export

DA.ds2 <- function(data, outcome, paired = NULL, p.adj = "fdr", ...){
  
  library(DESeq2)
  
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
  
  if(is.null(paired)){
    outcomedf <- data.frame(outcome = factor(outcome))
    row.names(outcomedf) <- colnames(count_table)
    x <- DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = outcomedf , design = ~ outcome)
  } else {
    outcomedf <- data.frame(outcome = factor(outcome),
                            paired = factor(paired))
    row.names(outcomedf) <- colnames(count_table)
    x <- DESeqDataSetFromMatrix(countData = as.data.frame(count_table), colData = outcomedf , design = ~ paired + outcome)
  }
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(x), 1, gm_mean)
  x = estimateSizeFactors(x, geoMeans = geoMeans)
  x <- DESeq(x, ...)
  res <- as.data.frame(results(x)@listData)
  colnames(res)[5] <- "pval"
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- results(x)@rownames
  res$Method <- "DESeq2"

  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  return(res)  
}

