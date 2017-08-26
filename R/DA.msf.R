#' MetagenomeSeq Feature model
#'
#' Implemented as in:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param outcome Factor. The outcome of interest. E.g. case and control
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param ... Additional arguments for the fitFeatureModel function

#' @export

DA.msf <- function(count_table, outcome, p.adj = "fdr", ...){

  library(metagenomeSeq, quietly = TRUE)
  count_table <- as.data.frame.matrix(count_table)
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  mod <- model.matrix(~outcome)
  mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgsdata,mod=mod,...)
  temp_table <- MRtable(mgsfit, number=nrow(count_table))
  temp_table <- temp_table[!is.na(row.names(temp_table)),]
  temp_table$Feature <- rownames(temp_table)
  colnames(temp_table)[7] <- "pval"
  temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
  temp_table$Method <- "MetagenomeSeq Feature"  

  return(temp_table)
}
