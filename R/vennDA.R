#' Plot Venn diagram from \code{allDA} object
#'
#' Plot a Venn (Euler) diagram of features found by different methods.
#' 
#' Require the eulerr package unless output is TRUE.
#' @param x (Required) Output from the \code{allDA} function
#' @param tests (Required) Character vector with tests to plot (E.g. \code{c("ttt","adx.t","wil")}, see \code{names(x$results)}). Default none
#' @param alpha Numeric. q-value threshold for significant features. Default 0.1
#' @param split If TRUE will split diagrams in positive and negative estimates if possible
#' @param output If TRUE will return a data.frame instead of a plot
#' @param pkg Use either "eulerr" package (default) or "venneuler" for drawing diagrams.
#' @param ... Additional arguments for plotting
#' @return If output TRUE then a data.frame with Features detected by the different methods
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(5)
#' mat <- matrix(rnbinom(500, size = 0.1, mu = 500), nrow = 50, ncol = 10)
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running allDA to compare methods
#' # This example uses 1 core (cores = 1). 
#' # Remove the cores argument to get it as high (and thereby fast) as possible.
#' res <- allDA(data = mat, predictor = pred, cores = 1)
#' 
#' # Plot venn diagram comparing significant features from znb and zpo
#' # znb and zpo only have significant features due to high false positive rates in this example
#' # split = TRUE splits the significant features in positive and negative estimates
#' vennDA(res, tests = c("znb","zpo"), split = TRUE)
#' @export
vennDA <- function(x, tests = NULL, alpha = 0.1, split = FALSE, output = FALSE, pkg = "eulerr", ...){

  # Load package
  if(pkg == "eulerr"){
    ok <- tryCatch({
      loadNamespace("eulerr")
      TRUE
    }, error=function(...) FALSE)
  } 
  if(pkg == "venneuler"){
    ok <- tryCatch({
      loadNamespace("venneuler")
      TRUE
    }, error=function(...) FALSE)
  }
  
  if (ok){
    # Check input
    if(!all(names(x) == c("raw","adj","est","details","results"))) stop("x is not an allDA object")
    
    plottests <- tests[tests %in% names(x[[2]])]  
    if(!all(tests %in% names(x[[2]]))){
      message(paste(tests[!tests %in% names(x[[2]])],collapse = ", ")," not found in the allDA object")
    }
    if(length(plottests) == 0) stop("Nothing to plot")
    
    # Which are significant
    featurelist <- list()
    for(i in seq_along(plottests)){
      sub <- x$adj[,c("Feature",plottests[i])]
      if(!plottests[i] %in% c("sam")) featurelist[[i]] <- sub[sub[,2] < alpha,"Feature"]
      if(plottests[i] %in% c("sam")) featurelist[[i]] <- sub[sub[,2] != "No","Feature"]
    }
    
    # Split in negative and positive significant
    if(split){
      featurelist.pos <- list()
      featurelist.neg <- list()
      for(i in seq_along(plottests)){
        
        subs <- x$est[,c(1,which(gsub("_.*","",colnames(x$est)) == plottests[i]))]
        
        if(plottests[i] == "bay"){
          sub.p <- subs[subs[,2] == levels(subs[,2])[1],"Feature"]
          sub.n <- subs[subs[,2] == levels(subs[,2])[2],"Feature"]
          featurelist.pos[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.p]
          featurelist.neg[[i]] <- featurelist[[i]][featurelist[[i]] %in% sub.n]
        }
        if(plottests[i] %in% c("abc","mva","sam","znb","zpo","poi","qpo","neb","lim","lli","lli2","vli","lia","lic","pea","spe","per","adx.t","adx.w","wil","ttt","ttr","ltt","ltt2","tta","ttc","ere","ere2","erq","erq2","ds2","ds2x","msf","zig","rai")){
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
        # If no estimate/logFC provided throw all significant in both positive and negative list
        if(!plottests[i] %in% c("abc","mva","sam","bay","znb","zpo","poi","qpo","neb","lim","lli","lli2","vli","lia","lic","pea","spe","per","adx.t","adx.w","wil","ttt","ttr","ltt","ltt2","tta","ttc","ere","ere2","erq","erq2","ds2","ds2x","msf","zig","rai")){
          featurelist.pos[[i]] <- featurelist[[i]]
          featurelist.neg[[i]] <- featurelist[[i]]
        }
      }
    } 
    
    # Collect significant features and make correct naming
    if(split){
      
      vennfeat.p <- do.call(c, featurelist.pos)
      vennfeat.n <- do.call(c, featurelist.neg)
      vennfeat <- c(vennfeat.p,vennfeat.n)
      if(length(vennfeat) == 0) stop("No significant features")
      naming.pos <- list()
      naming.neg <- list()
      for(i in seq_along(featurelist)){
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
      for(i in seq_along(featurelist)){
        naming[[i]] <- rep(plottests[i],length(featurelist[[i]]))
      }
      vennname <- do.call(c, naming)
      
    }
    
    # Make dataframe with significant features for each method
    venndf <- data.frame(vennfeat,vennname)
    
    # Remove the duplicate ones created earlier for methods without estimates/logFC
    if(split){
      for(i in seq_along(plottests)){
        if(!plottests[i] %in% c("abc","mva","sam","bay","znb","zpo","poi","qpo","neb","lim","lli","lli2","vli","lia","lic","pea","spe","per","adx.t","adx.w","wil","ttt","ltt","ltt2","tta","ttc","ere","ere2","erq","erq2","ds2","ds2x","msf","zig","rai")){
          venndf <- venndf[venndf$vennname != paste0(plottests[i],"_Negative"),]
          venndf$vennname <- as.character(venndf$vennname)
          venndf[venndf$vennname == paste0(plottests[i],"_Positive"),"vennname"] <- plottests[i]
        }
        if(plottests[i] %in% c("znb","zpo","poi","qpo","neb") & !plottests[i] %in% gsub("_.*","",colnames(x$est))){
          venndf <- venndf[venndf$vennname != paste0(plottests[i],"_Negative"),]
          venndf$vennname <- as.character(venndf$vennname)
          venndf[venndf$vennname == paste0(plottests[i],"_Positive"),"vennname"] <- plottests[i]
        }
      }
      
    }
    
    # Remove NAs
    venndf <- na.omit(venndf)
    
    # Return data.frame or plot
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
        plot(eulerr::euler(euler.list), quantities = TRUE, ...)
      }
    }
  } else {
    if(pkg == "eulerr") stop("eulerr package required")
    if(pkg == "venneuler") stop("venneuler package required")
  }

}
