#' Comparing differential abundance methods by FPR and AUC
#'
#' Calculating false positive rates and AUC for various differential abundance methods
#' @param otu_table Matrix or data.frame. OTU table with taxa as rows
#' @param predictor Factor. The outcome of interest. Should have two levels, e.g. case and control
#' @param R Integer. Number of times to run the models. Default 5
#' @param tests Character. Which tests to include. Default all
#' @param spikeMethod Character. Multiplicative ("mult") or additive ("add") spike-in. Default "mult"
#' @param effectSize Integer. The effect size for the spike-ins. Default 2
#' @param k Integer. Number of OTUs to spike in each tertile. k=5: 15 OTUs in total. Default 5
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param delta1 Numeric. The pseudocount for the Log t.test method. Default 1
#' @param delta2 Numeric. The pseudocount for the Log t.test2 method. Default 0.001
#' @param noOfIterations Integer. How many iterations should be run for the permutation test. Default 10000
#' @param margin Integer. The margin of when to stop iterating for non-significant OTUs for the permutation test. Default 50
#' @param testStat Function. The test statistic function for the permutation test (also in output of ttt, ltt, ltt2 and wil). Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param mc.samples Integer. Monte Carlo samples for ALDEx2. Default 64

#' @details Currently implemented methods:
#' \itemize{
#'  \item wil - Wilcoxon Rank Sum on relative abundances
#'  \item ttt - Welch t.test on relative abundances
#'  \item ltt - Welch t.test, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item ltt2 - Welch t.test, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item neb - Negative binomial GLM with log of library size as offset
#'  \item erq - EdgeR - Quasi likelihood
#'  \item ere - Edge R - Exact test
#'  \item msf - MetagenomeSeq feature model
#'  \item zig - MetagenomeSeq zero-inflated gaussian
#'  \item ds2 - DESeq2
#'  \item per - Permutation test with user defined test statistic
#'  \item bay - baySeq
#'  \item adx - ALDEx t-test and wilcoxon
#' }

#' @return A list of results:
#' \itemize{
#'  \item summary - A summary of the results
#'  \item table - FPR, AUC and spike detection rate for each run
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @import snow doSNOW foreach pROC
#' @importFrom parallel detectCores
#' @export

testDA <- function(otu_table, predictor, R = 5, tests = c("wil","ttt","ltt","ltt2","neb","erq","ere","msf","zig","ds2","per","bay","adx"), spikeMethod = "mult", effectSize = 2, k = 5, cores = (detectCores()-1), rng.seed = 123, delta1 = 1, delta2 = 0.001, noOfIterations = 10000, margin = 50, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, mc.samples = 64){

  if(sum(colSums(otu_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(otu_table) != length(predictor)) stop("Number of samples in OTU table does not match length of predictor")
  
  set.seed(rng.seed)
  
  final.results <- foreach(r = 1:R) %do% {
    
    library(parallel, quietly = TRUE)
    library(doSNOW, quietly = TRUE)
    library(foreach, quietly = TRUE)
    library(pROC, quietly = TRUE)
    
    message(paste0(r,". Run:"))
    
    # Shuffle predictor
    rand <- sample(predictor)
    
    # Remove OTUs not present in any samples
    otu_table <- otu_table[rowSums(otu_table) > 0,]
    
    # Spikein
    spiked <- spikein(otu_table, rand, spikeMethod, effectSize,  k)
    otu_table <- spiked[[1]]
    
    ### Run tests
    # Progress bar
    pb <- txtProgressBar(max = length(tests), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Start parallel
    if(cores == 1) {
      registerDoSEQ() 
    } else {
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
    
    results <- foreach(i = tests, .export = noquote(paste0("DA.",tests)), .options.snow = opts) %dopar% {
      
      res.sub <- switch(i,
                        wil = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,testStat)),
                        ttt = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,testStat)),
                        ltt = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,delta1,testStat)),
                        ltt2 = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,delta2,testStat)),
                        neb = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        erq = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        ere = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        msf = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        zig = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        ds2 = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        per = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,noOfIterations,rng.seed,margin,testStat)),
                        bay = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand)),
                        adx = do.call(get(noquote(paste0("DA.",i))),list(otu_table,rand,mc.samples)))
      
      
      return(res.sub)
      
    }
    if(cores != 1) stopCluster(cl)
    names(results) <- tests

    # Split ALDEx2 results in t.test and wilcoxon
    if("adx" %in% tests){
      adx.t <- as.data.frame(results["adx"])[,c(1:7,12)]
      adx.w <- as.data.frame(results["adx"])[,c(1:7,12)]
      colnames(adx.t) <- gsub("adx.","",colnames(adx.t))
      colnames(adx.w) <- colnames(adx.t)
      adx.t$pval <- as.numeric(as.data.frame(results["adx"])$adx.we.ep)
      adx.w$pval <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep)
      adx.t$Method <- "ALDEx2 t-test"
      adx.w$Method <- "ALDEx2 wilcoxon"
      results["adx"] <- NULL
      res.names <- names(results)
      results <- c(results,list(adx.t),list(adx.w))
      names(results) <- c(res.names,"adx.t","adx.w")
    }

    # Confusion matrix
    totalPos <- sapply(results,function(x) nrow(x[x$pval < 0.05,]))
    totalNeg <- sapply(results,function(x) nrow(x[x$pval >= 0.05,])) 
    trueNeg <- totalNeg  #if effectSize == 1
    truePos <- 0  #if effectSize == 1
    falseNeg <- 0 #if effectSize == 1
    if(effectSize != 1){
      truePos <- sapply(results, function(x) sum(x[x$pval < 0.05,"OTU"] %in% spiked[[2]]))
      falseNeg <- sapply(results, function(x) sum(x[x$pval >= 0.05,"OTU"] %in% spiked[[2]]))
    }
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR 
    fprs <- sapply(1:length(results), function(x) {
      if(totalPos[x] != 0){
        falsePos[x] / (totalPos[x] + totalNeg[x])
      } else {0}})
    
    
    # Spike detection rate
    sdrs <- sapply(1:length(results), function(x) truePos[x] / (k*3))
    
    # AUC
    aucs <- sapply(1:length(results), function(x) {
      if(effectSize != 1){
        test_roc <- NULL
        tryCatch(
          test_roc <- roc(as.numeric(results[[x]]$OTU %in% spiked[[2]]) ~ results[[x]]$pval, auc=TRUE, direction = ">"),
          error = function(e) NULL)
        if(!is.null(test_roc)){
          as.numeric(test_roc$auc) 
        }
      } else {
        0.5
      }
    })
    
    # Combine and return
    df.combined <- data.frame(Method = sapply(results, function(x) x$Method[1]),
                              FPR = fprs,
                              AUC = aucs,
                              Spike.detect.rate = sdrs,
                              Run = r)
    rownames(df.combined) <- NULL
    
    return(list(df.combined,results))
    
  }
  
  output.results <- do.call(rbind,lapply(final.results, function(x) x[[1]]))
  
  output.summary.fpr <- aggregate(FPR ~ Method, data = output.results, FUN = function(x) round(median(x),3))
  output.summary.auc <- aggregate(AUC ~ Method, data = output.results, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = output.results, FUN = function(x) round(median(x),3))
  output.summary <- merge(merge(output.summary.fpr,output.summary.auc, by = "Method"),output.summary.sdr, by = "Method")
    output.all.results <- lapply(final.results, function(x) x[[2]])
  
  return(list(summary = output.summary, table = output.results, results = output.all.results))
  
}




