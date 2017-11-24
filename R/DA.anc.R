#' ANCOM
#'
#' Implementation of ANCOM for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param allResults If TRUE will return raw results from the \code{ANCOM} function
#' @param ... Additional arguments for the \code{ANCOM} function
#' @export

DA.anc <- function(data, predictor, paired = NULL, allResults = FALSE, ...){
  
  suppressMessages(library(ancom.R, quietly = TRUE))

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
  } else {
    count_table <- data
  }
  
  # Ready data
  input <- as.data.frame(t(count_table))
  input$Group <- predictor
  if(!is.null(paired)) input$ID <- paired
  
  # Run test
  if(is.null(paired)){
    res <- ANCOM(input, ...) 
  } else {
    res <- ANCOM(input, repeated = TRUE, ...)
  }
  
  # Extract results
  df <- data.frame(Feature = rownames(count_table),
                   W = res[[1]],
                   Detected = factor("No",levels=c("No","Yes")))  
  df[df$Feature %in% res[[2]],"Detected"] <- "Yes"
  df$Method <- "ANCOM (anc)"
  
  if(class(data) == "phyloseq") df <- add.tax.DA(data, df)
  
  if(allResults){
    return(res)
  } else {
    return(df)
  }
  
}
