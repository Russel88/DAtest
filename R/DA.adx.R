#' Aldex t.test and wilcox
#' 
#' Implementation of \code{aldex} for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param ... Additional arguments for the \code{aldex} function
#' @export

DA.adx <- function(data, predictor, p.adj = "fdr", ...){
  
  suppressMessages(library(ALDEx2))
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
  } else {
    count_table <- data
  }

  # Run test
  x <- aldex(data.frame(count_table), as.character(predictor), ...)
  x <- x[,-c(9,11)]
  x$we.ep.adj <- p.adjust(x$we.ep, method = p.adj)
  x$wi.ep.adj <- p.adjust(x$wi.ep, method = p.adj)
  x$ordering <- NA
  x[!is.na(x$effect) & x$effect > 0,"ordering"] <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
  x[!is.na(x$effect) & x$effect < 0,"ordering"] <- paste0(levels(as.factor(predictor))[1],">",levels(as.factor(predictor))[2])
  x$Feature <- rownames(x)
  
  if(class(data) == "phyloseq") x <- add.tax.DA(data, x)
  
  return(x)
  
}




