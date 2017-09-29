#' SamSeq
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param fdr.output Passed to SAMseq. (Approximate) False Discovery Rate cutoff for output in significant genes table
#' @param allResults If TRUE will return raw results from the SAMseq function
#' @param ... Additional arguments for the SAMseq function
#' @import samr
#' @export
DA.sam <- function(data, predictor, paired = NULL, fdr.output = 0.05, allResults = FALSE, ...){

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
  } else {
    count_table <- data
  }

  if(is.numeric(predictor)){
    # Quantitative
    res <- samr::SAMseq(count_table, predictor, resp.type = "Quantitative", genenames = rownames(count_table), fdr.output = fdr.output, ...)
  } else {
    # Categorical
    predictor <- as.numeric(as.factor(predictor))
    
    if(length(levels(as.factor(predictor))) == 2){
      if(is.null(paired)){
        res <- samr::SAMseq(count_table, predictor, resp.type = "Two class unpaired", genenames = rownames(count_table), fdr.output = fdr.output)
      } else {
        predictor[predictor == 2] <- -1
        predictor <- as.numeric(as.factor(paired)) * predictor
        res <- samr::SAMseq(count_table, predictor, resp.type = "Two class paired", genenames = rownames(count_table), fdr.output = fdr.output, ...)
      }
    } else {
      res <- samr::SAMseq(count_table, predictor, resp.type = "Multiclass", genenames = rownames(count_table), fdr.output = fdr.output, ...)
    }
  }
  
  if(res$samr.obj$resp.type == "Multiclass"){
    df <- data.frame(Feature = rownames(count_table),
                     Score = res$samr.obj$tt,
                     Sig = factor("No",levels = c("No","Yes")))
    tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.up)[,1],"Sig"] <- "Yes",error = function(e) NULL)
    cont <- as.data.frame(res$samr.obj$stand.contrasts)
    colnames(cont) <- paste("Contrast",colnames(cont))
    df <- cbind(df,cont)
  } else {
    if(res$samr.obj$resp.type == "Quantitative"){
      df <- data.frame(Feature = rownames(count_table),
                       Score = res$samr.obj$tt,
                       Sig.up = factor("No",levels = c("No","Yes")),
                       Sig.lo = factor("No",levels = c("No","Yes")))
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.up)[,1],"Sig.up"] <- "Yes",error = function(e) NULL)
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.lo)[,1],"Sig.lo"] <- "Yes",error = function(e) NULL)
    } else {
      df <- data.frame(Feature = rownames(count_table),
                       Score = res$samr.obj$tt,
                       Fold.change = res$samr.obj$foldchange,
                       Sig.up = factor("No",levels = c("No","Yes")),
                       Sig.lo = factor("No",levels = c("No","Yes")))
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.up)[,1],"Sig.up"] <- "Yes",error = function(e) NULL)
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.lo)[,1],"Sig.lo"] <- "Yes",error = function(e) NULL)
    }
  }
  
  df$Method <- "SAMseq (sam)"
  
  if(class(data) == "phyloseq") df <- add.tax.DA(data, df)
  
  if(allResults){
    return(res)
  } else {
    return(df)
  }

}




