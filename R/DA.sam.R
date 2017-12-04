#' SAMSeq
#' 
#' SAMSeq implementation for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param fdr.output Passed to \code{SAMseq}. (Approximate) False Discovery Rate cutoff for output in significant genes table
#' @param allResults If TRUE will return raw results from the \code{SAMseq} function
#' @param ... Additional arguments for the \code{SAMseq} function
#' @export
DA.sam <- function(data, predictor, paired = NULL, fdr.output = 0.05, allResults = FALSE, ...){

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
  } else {
    count_table <- data
  }

  # Load package and make sure it unloads upon finishing to reset global error options
  suppressMessages(library(samr))
  on.exit(detach("package:samr", unload = TRUE))
  
  pred.lev <- levels(as.factor(predictor))
  
  # Run the test
  if(is.numeric(predictor)){
    # Quantitative
    res <- samr::SAMseq(count_table, predictor, resp.type = "Quantitative", genenames = rownames(count_table), fdr.output = fdr.output, ...)
  } else {
    # Categorical
    predictor <- as.numeric(as.factor(predictor))
    
    if(length(levels(as.factor(predictor))) == 2){
      if(is.null(paired)){
        res <- samr::SAMseq(count_table, predictor, resp.type = "Two class unpaired", genenames = rownames(count_table), fdr.output = fdr.output, ...)
      } else {
        predictor[predictor == 1] <- -1
        predictor[predictor == 2] <- 1
        predictor <- as.numeric(as.factor(paired)) * predictor
        res <- samr::SAMseq(count_table, predictor, resp.type = "Two class paired", genenames = rownames(count_table), fdr.output = fdr.output, ...)
      }
    } else {
      res <- samr::SAMseq(count_table, predictor, resp.type = "Multiclass", genenames = rownames(count_table), fdr.output = fdr.output, ...)
    }
  }
  
  # Collect results
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
      df$Sig <- "No"
      df[df$Sig.up == "Yes" | df$Sig.lo == "Yes","Sig"] <- "Yes"
    } else {
      df <- data.frame(Feature = rownames(count_table),
                       Score = res$samr.obj$tt,
                       Fold.change = res$samr.obj$foldchange,
                       log2FC = log2(res$samr.obj$foldchange),
                       Sig.up = factor("No",levels = c("No","Yes")),
                       Sig.lo = factor("No",levels = c("No","Yes")))
      df$ordering <- NA
      df[!is.na(df$log2FC) & df$log2FC > 0,"ordering"] <- paste0(pred.lev[2],">",pred.lev[1])
      df[!is.na(df$log2FC) & df$log2FC < 0,"ordering"] <- paste0(pred.lev[1],">",pred.lev[2])
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.up)[,1],"Sig.up"] <- "Yes",error = function(e) NULL)
      tryCatch(df[df$Feature %in% as.matrix(res$siggenes.table$genes.lo)[,1],"Sig.lo"] <- "Yes",error = function(e) NULL)
      df$Sig <- "No"
      df[df$Sig.up == "Yes" | df$Sig.lo == "Yes","Sig"] <- "Yes"
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




