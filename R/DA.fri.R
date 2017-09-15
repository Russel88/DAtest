#' Friedman Rank Sum test
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the friedman.test function
#' @param ... Additional arguments for the friedman.test function
#' @export

DA.fri <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  if(is.null(paired)) stop("Friedman test needs a paired argument")
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
  } else {
    count_table <- data
  }
  
  fri <- function(x){
    tryCatch(friedman.test(as.numeric(x), predictor, paired, ...), error = function(e){NA}) 
  }

  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  reslist <- apply(count.rel,1,fri)

  if(allResults){
    return(reslist)
  } else {
    res <- data.frame(statistic = sapply(reslist, function(x) x$statistic),
                      parameter = sapply(reslist, function(x) x$parameter),
                      pval = sapply(reslist, function(x) x$p.value))
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    
    res$Feature <- gsub(".Friedman.*","",rownames(res))
    res$Method <- "Friedman (fri)" 
    
    if(class(data) == "phyloseq"){
      if(!is.null(tax_table(data, errorIfNULL = FALSE))){
        tax <- tax_table(data)
        res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
        rownames(res) <- NULL
      } 
    }
    return(res)
  }
  
}
