#' ANCOM
#'
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param allResults If TRUE will return raw results from the ANCOM function
#' @param ... Additional arguments for the ANCOM function
#' @export

DA.anc <- function(data, predictor, paired = NULL, allResults = FALSE, ...){
  
  suppressMessages(library(ancom.R, quietly = TRUE))

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- as.factor(unlist(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
  } else {
    count_table <- data
  }
  
  input <- as.data.frame(t(count_table))
  input$Group <- predictor
  if(!is.null(paired)) input$ID <- paired
  
  if(is.null(paired)){
    res <- ANCOM(input, ...) 
  } else {
    res <- ANCOM(input, repeated = TRUE, ...)
  }
  
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
