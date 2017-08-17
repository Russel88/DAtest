#' Run many differential abundance methods
#'
#' Run many differential abundance tests at a time
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param predictor Factor. The outcome of interest. E.g. case and control
#' @param paired Factor. Subject ID for running paired analysis. Only for "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq" and "ds2"
#' @param relative Logical. Should abundances be made relative? Only has effect for "ttt", "wil", "per", "aov", "kru" and "lim". Default TRUE
#' @param tests Character. Which tests to include. Default all (See below for details)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for pvalue adjustment. Default "fdr"
#' @param delta1 Numeric. The pseudocount for the Log t.test/ANOVA/llm method. Default 1
#' @param delta2 Numeric. The pseudocount for the Log t.test2/ANOVA2/llm2 method. Default 0.001
#' @param noOfIterations Integer. How many iterations should be run for the permutation tests. Default 10000
#' @param margin Integer. The margin of when to stop iterating for non-significant Features for the permutation tests. Default 50
#' @param testStat Function. The test statistic function for the permutation test (also in output of ttt, ltt, ltt2 and wil). Should take two vectors as arguments. Default is a log fold change: log((mean(case abundances)+1)/(mean(control abundances)+1))
#' @param testStat.pair Function. The test statistic function for the paired permutation test (also in output of ttt, ltt, ltt2 and wil). Should take two vectors as arguments. Default is a log fold change: mean(log((case abundances+1)/(control abundances+1)))
#' @param mc.samples Integer. Monte Carlo samples for ALDEx2. Default 64
#' @param sig Numeric. Alpha used in ANCOM. Default 0.05
#' @param multcorr Integer. Correction used in ANCOM. Default 3 (no correction)
#' @param tau Numeric. Tuning parameter for ANCOM. Default 0.02
#' @param theta Numeric. Tuning parameter for ANCOM. Default 0.1
#' @param repeated Logical. Are there repeated measures? Only for ANCOM. Default FALSE
#' @param TMM.option 1 or 2. For "enn". Option "1" is for an approach using the mean of the	effective library sizes as a reference library size in TMM normalization; option "2" represents an approach to regenerating counts with a common dispersion. Default 1
#' @param log.lim Logical. For "lim". Should abundances be log transformed (before relative abundances if relative = TRUE). Default FALSE
#' @param delta.lim Pseudocount for "lim" log transformation. Default 1
#' @details Currently implemented methods:
#' \itemize{
#'  \item per - Permutation test with user defined test statistic
#'  \item bay - baySeq
#'  \item adx - ALDEx t-test and wilcoxon
#'  \item wil - Wilcoxon Rank Sum on relative abundances
#'  \item ttt - Welch t.test on relative abundances
#'  \item ltt - Welch t.test, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item ltt2 - Welch t.test, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item neb - Negative binomial GLM with log of library size as offset
#'  \item erq - EdgeR - Quasi likelihood
#'  \item ere - EdgeR - Exact test
#'  \item msf - MetagenomeSeq feature model
#'  \item zig - MetagenomeSeq zero-inflated gaussian
#'  \item ds2 - DESeq2
#'  \item enn - ENNB: Two-stage procedure from https://cals.arizona.edu/~anling/software.htm
#'  \item anc - ANCOM. This test does not output pvalues; for comparison with the other methods, detected Features are set to a pvalue of 0, all else are set to 1.
#'  \item lim - LIMMA. Moderated linear models based on emperical bayes
#'  \item kru - Kruskal-Wallis on relative abundances
#'  \item aov - ANOVA on relative abundances
#'  \item lao - ANOVA, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item lao2 - ANOVA, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item lrm - Linear regression on relative abundances
#'  \item llm - Linear regression, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item llm2 - Linear regression, but with relative abundances transformed with log(relative abundance + delta2)
#' }
#' Is it too slow? Remove "anc" from test argument.
#' 
#' "per" is also somewhat slow, but is usually one of the methods performing well.
#' @return A list of results:
#' \itemize{
#'  \item table - Summary of results
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @import snow doSNOW foreach pROC
#' @importFrom parallel detectCores
#' @export

allDA <- function(count_table, predictor, paired = NULL, relative = TRUE, tests = c("anc","per","bay","adx","enn","wil","ttt","ltt","ltt2","neb","erq","ere","msf","zig","ds2","lim","aov","lao","lao2","kru","lrm","llm","llm2"), cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", delta1 = 1, delta2 = 0.001, noOfIterations = 10000, margin = 50, testStat = function(case,control){log((mean(case)+1)/(mean(control)+1))}, testStat.pair = function(case,control){mean(log((case+1)/(control+1)))}, mc.samples = 64, sig = 0.05, multcorr = 3, tau = 0.02, theta = 0.1, repeated = FALSE, TMM.option = 1, log.lim = FALSE, delta.lim = 1){

  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(levels(as.factor(predictor))) < 2) stop("Predictor should have at least two levels")
  
  library(parallel, quietly = TRUE)
  library(doSNOW, quietly = TRUE)
  library(foreach, quietly = TRUE)
  library(pROC, quietly = TRUE)
  
  # Prune test argument
  if(!"baySeq" %in% rownames(installed.packages())) tests <- tests[tests != "bay"]
  if(!"ALDEx2" %in% rownames(installed.packages())) tests <- tests[tests != "adx"] 
  if(!"MASS" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb","enn")]
  if(!"lme4" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb")]
  if(!"edgeR" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ere","erq","enn")]
  if(!"metagenomeSeq" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("msf","zig")]
  if(!"DESeq2" %in% rownames(installed.packages())) tests <- tests[tests != "ds2"]
  if(!"ancom.R" %in% rownames(installed.packages())) tests <- tests[tests != "anc"]  
  if(!"glmnet" %in% rownames(installed.packages())) tests <- tests[tests != "enn"] 
  if(!"limma" %in% rownames(installed.packages())) tests <- tests[tests != "lim"] 
  
  if(!is.null(paired)){
    tests <- tests[!tests %in% c("bay","adx","anc","enn","ere","msf","zig","aov","lao","lao2","kru")]
  } 
  
  if(length(levels(as.factor(predictor))) > 2){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","aov","lao","lao2","kru","lrm","llm","llm2")]
  } else {
    tests <- tests[!tests %in% c("aov","lao","lao2","kru","lrm","llm","llm2")]
  }
  
  if(is.numeric(predictor)){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","lrm","llm","llm2")]
  }
  
  set.seed(rng.seed)
  message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  count_table <- count_table[rowSums(count_table) > 0,]
  
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
                      wil = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,testStat,testStat.pair,paired, p.adj, relative)),
                      ttt = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,testStat,testStat.pair,paired, p.adj, relative)),
                      ltt = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,delta1,testStat,testStat.pair,paired, p.adj)),
                      ltt2 = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,delta2,testStat,testStat.pair,paired, p.adj)),
                      neb = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj)),
                      erq = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj)),
                      ere = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj)),
                      msf = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj)),
                      zig = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj)),
                      ds2 = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj)),
                      per = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired,noOfIterations,rng.seed,margin,testStat,testStat.pair, p.adj, relative)),
                      bay = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj)),
                      adx = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,mc.samples, p.adj)),
                      enn = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,TMM.option,p.adj)),
                      anc = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,sig,multcorr, tau, theta, repeated)),
                      lim = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,p.adj,relative,paired,log.lim,delta.lim)),
                      kru = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj, relative)),
                      aov = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor, p.adj, relative)),
                      lao = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,delta1, p.adj)),
                      lao2 = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,delta2, p.adj)),
                      lrm = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj)),
                      llm = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj, delta1)),
                      llm2 = do.call(get(noquote(paste0("DA.",i))),list(count_table,predictor,paired, p.adj, delta2)))
    
    res.sub[is.na(res.sub$pval),"pval"] <- 1
    res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
    
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
    adx.t$pval.adj <- p.adjust(adx.t$pval, method = p.adj)
    adx.w$pval.adj <- p.adjust(adx.w$pval, method = p.adj)
    adx.t$Method <- "ALDEx2 t-test"
    adx.w$Method <- "ALDEx2 wilcox"
    results["adx"] <- NULL
    res.names <- names(results)
    results <- c(results,list(adx.t),list(adx.w))
    names(results) <- c(res.names,"adx.t","adx.w")
  }
  
  # Positives
  Pos.raw <- sapply(results,function(x) x[x$pval < 0.05,"Feature"])
  Pos.adj <- sapply(results,function(x) x[x$pval.adj < 0.05,"Feature"])

  features <- row.names(count_table)
  counts <- sapply(features, function(y) sum(unlist(sapply(Pos.raw, function(x) x %in% y)))) / length(tests)
  counts.adj <- sapply(features, function(y) sum(unlist(sapply(Pos.adj, function(x) x %in% y)))) / length(tests)
  
  # Combine and return
  df.combined <- data.frame(Feature = features,
                            Detect.rate.raw = counts,
                            Detect.rate.adj = counts.adj)
  rownames(df.combined) <- NULL
  
  return(list(table = df.combined,results = results))

}




