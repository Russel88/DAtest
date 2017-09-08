#' Wilcoxon Rank Sum and Signed Rank Test

#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param testStat Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param testStat.pair Function. Function for calculating fold change. Should take two vectors as arguments. Default is a log fold change: mean(log((case abundances+1)/(control abundances+1)))
#' @param ... Additional arguments for the wilcox.test function
#' @export

DA.wil <- function(data, predictor, paired = NULL, relative = TRUE, p.adj = "fdr", testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, ...){
 
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
    if(!is.null(paired)) paired <- suppressWarnings(as.matrix(sample_data(data)[,paired]))
  } else {
    count_table <- data
  }
  wil <- function(x){
    tryCatch(wilcox.test(x ~ predictor, ...)$p.value, error = function(e){NA}) 
  }

  if(!is.null(paired)){
    count_table <- count_table[,order(paired)]
    predictor <- predictor[order(paired)]
    testStat <- testStat.pair
    wil <- function(x){
      tryCatch(wilcox.test(x ~ predictor, paired = TRUE, ...)$p.value, error = function(e){NA}) 
    }
  }
  
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  res <- data.frame(pval = apply(count.rel,1,wil))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  # Teststat
  predictor.num <- as.numeric(as.factor(predictor))-1
  testfun <- function(x){
    case <- x[predictor.num==1]
    control <- x[predictor.num==0]
    testStat(case,control) 
  }
  
  res$FC <- apply(count.rel,1,testfun)
  
  res$Feature <- rownames(res)
  res$Method <- "Wilcox" 
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  return(res)
 
}
