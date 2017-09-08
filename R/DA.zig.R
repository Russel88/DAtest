#' MetageonomeSeq ZIG
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the fitZig function
#' @export

DA.zig <- function(data, predictor, p.adj = "fdr", ...){
  
  library(metagenomeSeq)
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
  } else {
    count_table <- data
  }
  
  count_table <- as.data.frame.matrix(count_table)
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  mod <- model.matrix(~predictor)
  mgsfit <- fitZig(obj=mgsdata,mod=mod, ...)
  temp_table <- MRtable(mgsfit, number=nrow(count_table))
  temp_table <- temp_table[!is.na(row.names(temp_table)),]
  # Pvalue have different naming depending on package version
  #if("pvalues" %in% names(temp_table)) res.zig <- data.frame(pval = temp_table$pvalues)
  #if("pValue" %in% names(temp_table)) res.zig <- data.frame(pval = temp_table$pValue)
  colnames(temp_table)[8] <- "pval"
  temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
  temp_table$Feature <- rownames(temp_table)
  temp_table$Method <- "MetagenomeSeq ZIG"
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      temp_table <- merge(temp_table, tax, by.x = "Feature", by.y = "row.names")
      rownames(temp_table) <- NULL
    } 
  }
  
  return(temp_table)
}


