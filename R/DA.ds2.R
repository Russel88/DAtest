#' DESeq2
#'
#' From:
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0208-8
#' Manual geometric means calculated to avoid errors, see https://github.com/joey711/phyloseq/issues/387

#' @export

DA.ds2 <- function(count_table, outcome, paired = NULL, p.adj){
  
  library(DESeq2, quietly = TRUE)
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
  x <- DESeq(x)
  res <- as.data.frame(results(x)@listData)
  colnames(res)[5] <- "pval"
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  res$Feature <- results(x)@rownames
  res$Method <- "DESeq2"

  return(res)  
}

