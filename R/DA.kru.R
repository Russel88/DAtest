#' Kruskal-Wallis test
#' 
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Factor, OR if data is a phyloseq object the name of the variable in sample_data in quotation
#' @param relative Logical. Should count_table be normalized to relative abundances. Default TRUE
#' @param p.adj Character. P-value adjustment. Default "fdr". See p.adjust for details
#' @param allResults If TRUE will return raw results from the kruskal.test function
#' @param ... Additional arguments for the kruskal.test function
#' @export

DA.kru <- function(data, predictor, relative = TRUE, p.adj = "fdr", allResults = FALSE, ...){
 
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1) stop("When data is a phyloseq object predictor should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,predictor])))
  } else {
    count_table <- data
  }
  
  kru <- function(x){
    tryCatch(kruskal.test(as.numeric(x) ~ predictor, ...)$p.value, error = function(e){NA}) 
  }

  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }
  
  res <- data.frame(pval = apply(count.rel,1,kru))
  res$pval.adj <- p.adjust(res$pval, method = p.adj)
  
  res$Feature <- rownames(res)
  res$Method <- "Kruskal-Wallis (kru)" 
  
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      tax <- tax_table(data)
      res <- merge(res, tax, by.x = "Feature", by.y = "row.names")
      rownames(res) <- NULL
    } 
  }
  
  if(allResults){
    if(is.null(paired)){
      kru <- function(x){
        tryCatch(kruskal.test(x ~ predictor, ...), error = function(e){NA}) 
      }
    } else {
      kru <- function(x){
        tryCatch(kruskal.test(x ~ predictor, paired = TRUE, ...), error = function(e){NA}) 
      }
    }
    return(apply(count.rel,1,kru))
  } else {
    return(res)
  }
 
}
