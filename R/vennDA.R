#' Plot Venn diagram from allDA object
#'
#' Plot a Venn (Euler) diagram of features found by different methods.
#' 
#' Require the eulerr package unless output is TRUE.
#' @param x Output from the allDA function
#' @param tests Character vector with tests to plot (E.g. c("ttt","adx.t","wil"), see names(x$results)). Default none
#' @param alpha Numeric. q-value threshold for significant features. Default 0.05
#' @param split If TRUE will split diagrams in positive and negative estimates if possible
#' @param output If TRUE will return a data.frame instead of a plot
#' @param pkg Use either "eulerr" package (default) or "venneuler" for drawing diagrams.
#' @param ... Additional arguments for plotting
#' @return If output TRUE then a data.frame with Features detected by the different methods
#' @export
vennDA <- function(x, tests = NULL, alpha = 0.05, split = FALSE, output = FALSE, pkg = "eulerr", ...){

  if(pkg == "eulerr") library(eulerr)
  if(pkg == "venneuler") library(venneuler)
  
  if(!all(names(x) == c("raw","adj","est","results"))) stop("x is not an allDA object")
  
  plottests <- tests[tests %in% names(x[[2]])]  
  if(!all(tests %in% names(x[[2]]))){
    message(paste(tests[!tests %in% names(x[[2]])],collapse = ", ")," not found in the allDA object")
  }
  if(length(plottests) == 0) stop("Nothing to plot")
  
  featurelist <- list()
  for(i in seq_along(plottests)){
    sub <- x$adj[,c("Feature",plottests[i])]
    featurelist[[i]] <- sub[sub[,2] < alpha,"Feature"]
    if(plottests[i] == "sam") featurelist[[i]] <- sub[sub != "No","Feature"]
    if(plottests[i] == "anc") featurelist[[i]] <- sub[sub == "Yes","Feature"]
  }
  
  if(split){
    featurelist.pos <- list()
    featurelist.neg <- list()
    for(i in seq_along(plottests)){
      
      subs <- x$est[,c(1,which(gsub("_.*","",colnames(x$est)) == plottests[i]))]

      if(plottests[i] == "sam"){
        sub.p <- subs[subs[,2] > 1,"Feature"]
        sub.n <- subs[subs[,2] < 1,"Feature"]
        featurelist.pos[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.p]
        featurelist.neg[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.n]
      }
      if(plottests[i] == "bay"){
        sub.p <- subs[subs[,2] == levels(subs[,2])[1],"Feature"]
        sub.n <- subs[subs[,2] == levels(subs[,2])[2],"Feature"]
        featurelist.pos[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.p]
        featurelist.neg[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.n]
      }
      if(plottests[i] %in% c("znb","zpo","poi","qpo","neb","lrm","llm","llm2","lim","lli","lli2","vli","pea","spe","per","adx.t","adx.w","wil","ttt","ltt","ltt2","ere","ere2","erq","erq2","ds2","msf","zig")){
        if(is.null(ncol(subs))){
          featurelist.pos[[i]] <- featurelist[[i]]
          featurelist.neg[[i]] <- featurelist[[i]]
        } else {
          sub.p <- subs[subs[,2] > 0,"Feature"]
          sub.n <- subs[subs[,2] < 0,"Feature"]
          featurelist.pos[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.p]
          featurelist.neg[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.n]
        }
      }
      if(!plottests[i] %in% c("sam","bay","znb","zpo","poi","qpo","neb","lrm","llm","llm2","lim","lli","lli2","vli","pea","spe","per","adx.t","adx.w","wil","ttt","ltt","ltt2","ere","ere2","erq","erq2","ds2","msf","zig")){
        featurelist.pos[[i]] <- featurelist[[i]]
        featurelist.neg[[i]] <- featurelist[[i]]
      }
    }
  } 

  if(split){
    vennfeat.p <- do.call(c, featurelist.pos)
    vennfeat.n <- do.call(c, featurelist.neg)
    vennfeat <- c(vennfeat.p,vennfeat.n)
    if(length(vennfeat) == 0) stop("No significant features")
    naming.pos <- list()
    naming.neg <- list()
    for(i in 1:length(featurelist)){
      if(plottests[i] == "bay"){
        naming.pos[[i]] <- rep(paste0(plottests[i],"_",levels(x$results$bay$ordering)[1]),length(featurelist.pos[[i]]))
        naming.neg[[i]] <- rep(paste0(plottests[i],"_",levels(x$results$bay$ordering)[2]),length(featurelist.neg[[i]]))
      } else {
        naming.pos[[i]] <- rep(paste0(plottests[i],"_Positive"),length(featurelist.pos[[i]]))
        naming.neg[[i]] <- rep(paste0(plottests[i],"_Negative"),length(featurelist.neg[[i]])) 
      }
    }
    vennname.pos <- do.call(c, naming.pos)
    vennname.neg <- do.call(c, naming.neg)
    vennname <- c(vennname.pos,vennname.neg)

  } else {
    vennfeat <- do.call(c, featurelist)
    if(length(vennfeat) == 0) stop("No significant features")
    naming <- list()
    for(i in 1:length(featurelist)){
      naming[[i]] <- rep(plottests[i],length(featurelist[[i]]))
    }
    vennname <- do.call(c, naming)
  }
  
  venndf <- data.frame(vennfeat,vennname)

  for(i in seq_along(plottests)){
    if(!plottests[i] %in% c("sam","bay","znb","zpo","poi","qpo","neb","lrm","llm","llm2","lim","lli","lli2","vli","pea","spe","per","adx.t","adx.w","wil","ttt","ltt","ltt2","ere","ere2","erq","erq2","ds2","msf","zig")){
      venndf <- venndf[venndf$vennname != paste0(plottests[i],"_Negative"),]
      venndf$vennname <- as.character(venndf$vennname)
      venndf[venndf$vennname == paste0(plottests[i],"_Positive"),"vennname"] <- plottests[i]
    }
    if(plottests[i] %in% c("znb","zpo","poi","qpo","neb","lrm","llm","llm2") & !plottests[i] %in% gsub("_.*","",colnames(x$est))){
      venndf <- venndf[venndf$vennname != paste0(plottests[i],"_Negative"),]
      venndf$vennname <- as.character(venndf$vennname)
      venndf[venndf$vennname == paste0(plottests[i],"_Positive"),"vennname"] <- plottests[i]
    }
  }
  
  if(output){
    colnames(venndf) <- c("Feature","Method")
    return(venndf)
  } else {
    if(pkg == "venneuler"){
      venndia <- venneuler::venneuler(venndf)
      plot(venndia, ...)
    }
    if(pkg == "eulerr"){
      euler.list <- list()
      for(i in seq_along(unique(venndf$vennname))){
        euler.list[[i]] <- venndf[venndf$vennname == unique(venndf$vennname)[i],1]
      }
      names(euler.list) <- unique(venndf$vennname)
      plot(euler(euler.list), counts = TRUE, ...)
    }
  }
}
