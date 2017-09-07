#' Comparing differential abundance methods by FPR and AUC
#'
#' Calculating false positive rates and AUC (Area Under the receiver operator Curve) for various differential abundance methods
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param outcome The outcome of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the outcome is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli" and "lli2"
#' @param R Integer. Number of times to run the tests. Default 10
#' @param tests Character. Which tests to include. Default all (See below for details)
#' @param relative Logical. Should abundances be made relative? Only has effect for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm" and "spe". Default TRUE
#' @param effectSize Integer. The effect size for the spike-ins. Default 2
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). k=c(5,10,15): 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default c(5,5,5)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
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
#'  \item lim - LIMMA. Moderated linear models based on emperical bayes
#'  \item lli - LIMMA, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item lli2 - LIMMA, but with relative abundances transformed with log(relative abundance + delta2)
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
#' "neb" can be slow if there is a paired argument.
#' 
#' "per" is also somewhat slow, but is usually one of the methods performing well with large sample sizes.
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
#' @return An object of class DA, which contains a list of results:
#' \itemize{
#'  \item table - FPR, AUC and spike detection rate for each run
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @import snow doSNOW foreach
#' @importFrom parallel detectCores
#' @importFrom pROC roc
#' @export

testDA <- function(data, outcome, paired = NULL, R = 10, tests = c("neb","rai","per","bay","adx","wil","ttt","ltt","ltt2","erq","ere","msf","zig","ds2","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe"), relative = TRUE, effectSize = 2, k = c(5,5,5), cores = (detectCores()-1), rng.seed = 123, args = list(), verbose = FALSE){

  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(outcome) > 1 | length(paired) > 1) stop("When data is a phyloseq object outcome and paired should only contain the name of the variables in sample_data")
    if(!outcome %in% sample_variables(data)) stop(paste(outcome,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    outcome <- suppressWarnings(as.matrix(sample_data(data)[,outcome]))
    if(!is.null(paired)) paired <- suppressWarnings(as.matrix(sample_data(data)[,paired]))
  } else {
    count_table <- data
  }
  
  # Coerce data
  if(!is.null(paired)) paired <- as.factor(paired)
  count_table <- as.matrix(count_table)
  
  # Checks
  if(min(count_table) < 0) stop("Count_table contains negative values!")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(outcome)) stop("Number of samples in count_table does not match length of outcome")
  if(length(levels(as.factor(outcome))) < 2) stop("outcome should have at least two levels")
  
  # Prune tests argument
  tests <- prune.tests.DA(tests, outcome, paired, relative)
  
  if(verbose){
    message(paste("Tests are run in the following order:"))
    print(as.data.frame(tests))
  } 
  
  set.seed(rng.seed)
  if(verbose) message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  if(verbose) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Numeric outcome
  if(is.numeric(outcome)){
    num.pred <- TRUE
    print("outcome is assumed to be numeric")
  } else {
    num.pred <- FALSE
  }
  
  final.results <- foreach::foreach(r = 1:R) %do% {

    message(paste0(r,". Run:"))
    
    # Shuffle outcome
    if(is.null(paired)){
      rand <- sample(outcome)
    } else {
      rand <- unsplit(lapply(split(outcome,paired), sample), paired)
    }
    
    # Spikein
    spiked <- spikein(count_table, rand, effectSize,  k, num.pred, relative)
    count_table <- spiked[[1]]
    
    # Run tests
    results <- run.tests.DA(count_table, rand, paired, tests, relative, args, cores)
    
    # Insert spiked column
    newnames <- names(results)
    results <- foreach(rsp = 1:length(results)) %do% {
      temp <- results[[rsp]]
      temp$Spiked <- "No"
      temp[temp$Feature %in% spiked[[2]],"Spiked"] <- "Yes"
      return(temp)
    }
    names(results) <- newnames
    
    # Confusion matrix
    totalPos <- sapply(results,function(x) nrow(x[x$pval < 0.05,]))
    totalNeg <- sapply(results,function(x) nrow(x[x$pval >= 0.05,])) 
    trueNeg <- totalNeg  #if effectSize == 1
    truePos <- 0  #if effectSize == 1
    falseNeg <- 0 #if effectSize == 1
    if(effectSize != 1){
      truePos <- sapply(results, function(x) sum(x[x$pval < 0.05,"Feature"] %in% spiked[[2]]))
      falseNeg <- sapply(results, function(x) sum(x[x$pval >= 0.05,"Feature"] %in% spiked[[2]]))
    }
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR 
    fprs <- sapply(1:length(results), function(x) {
      if(totalPos[x] != 0){
        falsePos[x] / (totalPos[x] + totalNeg[x])
      } else {0}})
    
    
    # Spike detection rate
    sdrs <- sapply(1:length(results), function(x) truePos[x] / sum(k))
    
    # AUC
    aucs <- sapply(1:length(results), function(x) {
      if(effectSize != 1){
        test_roc <- NULL
        tryCatch(
          test_roc <- pROC::roc(as.numeric(results[[x]]$Feature %in% spiked[[2]]) ~ results[[x]]$pval, auc=TRUE, direction = ">"),
          error = function(e) NULL)
        if(!is.null(test_roc)){
          as.numeric(test_roc$auc) 
        } else {
          0.5
        }
      } else {
        0.5
      }
    })

    # Combine and return
    df.combined <- data.frame(Method = sapply(results, function(x) x$Method[1]),
                              AUC = aucs,
                              FPR = fprs,
                              Spike.detect.rate = sdrs,
                              Run = r)
    rownames(df.combined) <- NULL
    
    return(list(df.combined,results))
    
  }
  
  output.results <- do.call(rbind,lapply(final.results, function(x) x[[1]]))
  output.all.results <- lapply(final.results, function(x) x[[2]])
  
  out <- list(table = output.results, results = output.all.results)
  class(out) <- "DA"
  return(out)
}




