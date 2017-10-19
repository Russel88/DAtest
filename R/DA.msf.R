#' MetagenomeSeq Feature model
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the fitFeatureModel function
#' @param ... Additional arguments for the fitFeatureModel function
#' @export

DA.msf <- function(data, predictor, p.adj = "fdr", allResults = FALSE, ...){

  suppressMessages(library(metagenomeSeq))
  
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
  mgsfit <- fitFeatureModel(obj=mgsdata,mod=mod,...)
  temp_table <- MRtable(mgsfit, number=nrow(count_table))
  temp_table <- temp_table[!is.na(row.names(temp_table)),]
  temp_table$Feature <- rownames(temp_table)
  colnames(temp_table)[7] <- "pval"
  temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
  temp_table$Method <- "MgSeq Feature (msf)"  

  if(class(data) == "phyloseq") temp_table <- add.tax.DA(data, temp_table)
  
  if(allResults) return(mgsfit) else return(temp_table)
}
