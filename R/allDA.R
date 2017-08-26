#' Run many differential abundance methods
#'
#' Run many differential abundance tests at a time
#' @param count_table Matrix or data.frame. Table with taxa/genes/proteins as rows and samples as columns
#' @param predictor Factor or Numeric. The outcome of interest. E.g. case and control. If the predictor has more than two levels, only the 2. level will be spiked. If the predictor is numeric it will be treated as such in the analyses
#' @param paired Factor. Subject ID for running paired analysis. Only for "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli" and "lli2"
#' @param tests Character. Which tests to include. Default all (See below for details)
#' @param relative Logical. Should abundances be made relative? Only has effect for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm" and "spe". Default TRUE
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for pvalue adjustment. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param verbose Logical. Print information during run
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
#'  \item anc - ANCOM. This test does not output pvalues; for comparison with the other methods, detected Features are set to a pvalue of 0, all else are set to 1.
#'  \item lim - LIMMA. Moderated linear models based on emperical bayes
#'  \item kru - Kruskal-Wallis on relative abundances
#'  \item aov - ANOVA on relative abundances
#'  \item lao - ANOVA, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item lao2 - ANOVA, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item lrm - Linear regression on relative abundances
#'  \item llm - Linear regression, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item llm2 - Linear regression, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item rai - RAIDA
#'  \item spe - Spearman correlation
#' }
#' 
#' Additional arguments can be passed to the internal functions with the "args" argument. 
#' It should be structured as a list with elements named by the tests: 
#' E.g. passing to the DA.per function that it should only run 1000 iterations: args = list(per=list(noOfIterations=1000)).
#' Include that the log t.test should use a pseudocount of 0.1: args = list(per=list(noOfIterations=1000), ltt=list(delta=0.1)). 
#' Additional arguments are simply seperated by commas.
#' 
#' Below is an overview of which functions get the arguments that are passed to a specific test
#' \itemize{
#'  \item per - Passed to DA.per
#'  \item bay - Passed to getPriors and getLikelihoods
#'  \item adx - Passed to aldex
#'  \item wil - Passed to wilcox.test and DA.wil
#'  \item ttt - Passed to t.test and DA.ttt
#'  \item ltt - Passed to t.test and DA.ltt
#'  \item ltt2 - Passed to t.test and DA.ltt2
#'  \item neb - Passed to glm.nb and glmer.nb
#'  \item erq - Passed to exactTest
#'  \item ere - Passed to glmQLFit
#'  \item msf - Passed to fitFeatureModel
#'  \item zig - Passed to fitZig
#'  \item ds2 - Passed to DESeq
#'  \item anc - Passed to ANCOM
#'  \item lim - Passed to eBayes
#'  \item lli - Passed to eBayes
#'  \item lli2 - Passed to eBayes
#'  \item kru - Passed to kruskal.test
#'  \item aov - Passed to aov
#'  \item lao - Passed to aov
#'  \item lao2 - Passed to aov
#'  \item lrm - Passed to lm and lme
#'  \item llm - Passed to lm and lme
#'  \item llm2 - Passed to lm and lme
#'  \item rai - Passed to raida
#'  \item spe - Passed to cor.test
#' }
#' @return A list of results:
#' \itemize{
#'  \item table - Summary of results
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @import snow doSNOW foreach pROC
#' @importFrom parallel detectCores
#' @export

allDA <- function(count_table, predictor, paired = NULL, tests = c("spe","anc","per","bay","adx","wil","ttt","ltt","ltt2","neb","erq","ere","msf","zig","ds2","lim","aov","lao","lao2","kru","lrm","llm","llm2","rai"), relative = TRUE, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list(), verbose = TRUE){

  # Checks
  if(min(count_table) < 0 & is.numeric(predictor)) stop("Numeric predictor and negative values in count_table is currently not supported")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(levels(as.factor(predictor))) < 2) stop("Predictor should have at least two levels")
  
  library(parallel, quietly = TRUE)
  library(doSNOW, quietly = TRUE)
  library(foreach, quietly = TRUE)
  library(pROC, quietly = TRUE)
  
  # Prune test argument if packages are not installed
  if(!"baySeq" %in% rownames(installed.packages())) tests <- tests[tests != "bay"]
  if(!"ALDEx2" %in% rownames(installed.packages())) tests <- tests[tests != "adx"] 
  if(!"MASS" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb")]
  if(!"lme4" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb")]
  if(!"edgeR" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ere","erq")]
  if(!"metagenomeSeq" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("msf","zig")]
  if(!"DESeq2" %in% rownames(installed.packages())) tests <- tests[tests != "ds2"]
  if(!"ancom.R" %in% rownames(installed.packages())) tests <- tests[tests != "anc"]  
  if(!"limma" %in% rownames(installed.packages())) tests <- tests[tests != "lim"]
  if(!"RAIDA" %in% rownames(installed.packages())) tests <- tests[tests != "rai"]
  
  # Excluded tests that do not work with a paired argument
  if(!is.null(paired)){
    tests <- tests[!tests %in% c("bay","adx","anc","ere","msf","zig","aov","lao","lao2","kru","rai","spe")]
  } 
  
  # Only include some tests if there are more than two levels in predictor
  if(length(levels(as.factor(predictor))) > 2){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe")]
  } else {
    # Excluded tests if levels in predictor is exactly 2
    tests <- tests[!tests %in% c("aov","lao","lao2","kru","lrm","llm","llm2","spe")]
  }
  
  # Only include specific tests if predictor is numeric
  if(is.numeric(predictor)){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","lrm","llm","llm2","spe")]
  } else {
    # Exclude if not numeric
    tests <- tests[!tests %in% c("spe")]
  }
  
  # Exclude if relative is false
  if(relative == FALSE){
    tests <- tests[!tests %in% c("ltt2","neb","erq","ere","msf","zig","bay","ds2","adx","anc","lli2","lao2","llm2","rai")]
  }
  
  if(verbose){
    message(paste("Tests are run in the following order:"))
    print(as.data.frame(tests))
  } 
  
  set.seed(rng.seed)
  if(verbose) message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  if(verbose) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
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
    cl <- makeCluster(cores, outfile = "")
    registerDoSNOW(cl)
  }
  
  results <- foreach(i = tests, .export = noquote(paste0("DA.",tests)), .options.snow = opts) %dopar% {
    
    # Extract test arguments
    if(!all(names(args) %in% tests)) stop("One or more names in list with additional arguments does not match names of tests")
    for(j in seq_along(args)){
      assign(paste0(names(args)[j],".args"),args[[j]],pos=1)
    }
    test.args <- paste0(tests,".args")
    test.boo <- lapply(test.args,exists)
    for(l in seq_along(test.args)){
      if(test.boo[l] == FALSE) assign(test.args[l], list(),pos=1)
    }
    
    res.sub <- switch(i,
                      wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, relative, p.adj),wil.args)),
                      ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, relative, p.adj),ttt.args)),
                      ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired,relative, p.adj),ltt.args)),
                      ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),ltt2.args)),
                      neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),neb.args)),
                      erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),erq.args)),
                      ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, p.adj),ere.args)),
                      msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, p.adj),msf.args)),
                      zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, p.adj),zig.args)),
                      ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),ds2.args)),
                      per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, relative, p.adj),per.args)),
                      bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),bay.args)),
                      adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand),adx.args)),
                      anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand),anc.args)),
                      lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired,relative,p.adj),lim.args)),
                      lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired,relative,p.adj),lli.args)),
                      lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired,p.adj),lli2.args)),
                      kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, relative, p.adj),kru.args)),
                      aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, relative, p.adj),aov.args)),
                      lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,relative, p.adj),lao.args)),
                      lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand, p.adj),lao2.args)),
                      lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, relative, p.adj),lrm.args)),
                      llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired,relative, p.adj),llm.args)),
                      llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,paired, p.adj),llm2.args)),
                      rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,p.adj),rai.args)),
                      spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,rand,relative,p.adj),spe.args)))
    
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
  counts <- sapply(features, function(y) sum(unlist(sapply(Pos.raw, function(x) x %in% y)))) / length(results)
  counts.adj <- sapply(features, function(y) sum(unlist(sapply(Pos.adj, function(x) x %in% y)))) / length(results)
  
  # Combine and return
  df.combined <- data.frame(Feature = features,
                            Detect.rate.raw = counts,
                            Detect.rate.adj = counts.adj)
  rownames(df.combined) <- NULL
  
  return(list(table = df.combined,results = results))

}




