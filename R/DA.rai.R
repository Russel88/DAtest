#' RAIDA
#' 
#' Implementation of \code{raida} for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param ... Additional arguments for the \code{raida} function
#' @export

DA.rai <- function(data, predictor, p.adj = "fdr", ...){
  
  suppressMessages(library(RAIDA))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
  } else {
    count_table <- data
  }

  # Order count_table
  count_table.o <- as.data.frame(count_table[,order(predictor)])
  
  # Run test and collect results
  res <- raida(count_table.o, n.lib = as.numeric(table(predictor)), mtcm = p.adj, ...)
  res$log2FC <- log2(exp(res$mean2)/exp(res$mean1))
  res$ordering <- NA
  res[!is.na(res$log2FC) & res$log2FC > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
  res[!is.na(res$log2FC) & res$log2FC < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
  res$Feature <- rownames(res)
  colnames(res)[1] <- "pval"
  colnames(res)[2] <- "pval.adj"
  res$Method <- "RAIDA (rai)" 
  
  if(class(data) == "phyloseq") res <- add.tax.DA(data, res)
  
  return(res)

}




