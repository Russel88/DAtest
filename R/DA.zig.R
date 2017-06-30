#' MetageonomeSeq ZIG
#'
#' From
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8

#' @export

DA.zig <- function(count_table, outcome, p.adj){
  
  library(metagenomeSeq, quietly = TRUE)
  
  count_table <- as.data.frame.matrix(count_table)
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  mod <- model.matrix(~outcome)
  mgsfit <- fitZig(obj=mgsdata,mod=mod)
  temp_table <- MRtable(mgsfit, number=nrow(count_table))
  temp_table <- temp_table[!is.na(row.names(temp_table)),]
  # Pvalue have different naming depending on package version
  #if("pvalues" %in% names(temp_table)) res.zig <- data.frame(pval = temp_table$pvalues)
  #if("pValue" %in% names(temp_table)) res.zig <- data.frame(pval = temp_table$pValue)
  colnames(temp_table)[8] <- "pval"
  temp_table$pval.adj <- p.adjust(temp_table$pval, method = p.adj)
  temp_table$Feature <- rownames(temp_table)
  temp_table$Method <- "MetagenomeSeq ZIG"
  return(temp_table)
}


