#' MetagenomeSeq Feature model
#'
#' From:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.msf <- function(otu_table, outcome, p.adj){

  library(metagenomeSeq, quietly = TRUE)
  otu_table <- as.data.frame.matrix(otu_table)
  mgsdata <- newMRexperiment(counts = otu_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  mod <- model.matrix(~outcome)
  mgsfit <- metagenomeSeq::fitFeatureModel(obj=mgsdata,mod=mod)
  temp_table <- MRtable(mgsfit, number=nrow(otu_table))
  temp_table <- temp_table[!is.na(row.names(temp_table)),]
  temp_table$OTU <- rownames(temp_table)
  colnames(temp_table)[7] <- "pval"
  temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
  temp_table$Method <- "MetagenomeSeq Feature"  

  return(temp_table)
}
