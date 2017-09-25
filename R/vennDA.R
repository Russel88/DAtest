#' Plot Venn diagram from allDA object
#'
#' Plot a Venn (Euler) diagram of features found by different methods
#' @param x Output from the allDA function
#' @param tests Character vector with tests to plot (E.g. c("ttt","adx.t","wil"), see names(x$results)). Default none
#' @param alpha Numeric. q-value threshold for significant features. Default 0.05
#' @param ... Additional arguments for plotting
#' @return Nothing
#' @import venneuler
#' @export
vennDA <- function(x, tests = NULL, alpha = 0.05, ...){
  
  if(!all(names(x) == c("table","results"))) stop("x is not an allDA object")
  
  plottests <- tests[tests %in% names(x[[2]])]  
  if(!all(tests %in% names(x[[2]]))){
    message(paste(tests[!tests %in% names(x[[2]])],collapse = ",")," not found in the allDA object")
  }
  if(length(plottests) == 0) stop("Nothing to plot")
  
  featurelist <- list()
  for(i in 1:length(plottests)){
    sub <- x[[2]][[plottests[i]]]
    featurelist[[i]] <- sub[sub$pval.adj < alpha & !is.na(sub$pval.adj),"Feature"]
  }
  
  vennfeat <- do.call(c, featurelist)
  if(length(vennfeat) == 0) stop("No significant features")
  naming <- list()
  for(i in 1:length(featurelist)){
    naming[[i]] <- rep(plottests[i],length(featurelist[[i]]))
  }
  vennname <- do.call(c, naming)
  
  venndf <- data.frame(vennfeat,vennname)
  venndia <- venneuler(venndf)
  plot(venndia, ...)
  
}
