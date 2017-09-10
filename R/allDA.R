#' Run many differential abundance methods
#'
#' Run many differential abundance tests at a time
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param tests Character. Which tests to include. Default all (See below for details)
#' @param relative Logical. Should abundances be made relative? Only has effect for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea". Default TRUE
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for pvalue adjustment. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
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
#'  \item kru - Kruskal-Wallis on relative abundances
#'  \item aov - ANOVA on relative abundances
#'  \item lao - ANOVA, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item lao2 - ANOVA, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item lrm - Linear regression on relative abundances
#'  \item llm - Linear regression, but reads are first transformed with log(abundance + delta1) then turned into relative abundances
#'  \item llm2 - Linear regression, but with relative abundances transformed with log(relative abundance + delta2)
#'  \item rai - RAIDA
#'  \item spe - Spearman correlation
#'  \item pea - Pearson correlation
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
#'  \item lim - Passed to eBayes
#'  \item lli - Passed to eBayes and DA.lli
#'  \item lli2 - Passed to eBayes and DA.lli
#'  \item kru - Passed to kruskal.test
#'  \item aov - Passed to aov
#'  \item lao - Passed to aov and DA.lao
#'  \item lao2 - Passed to aov and DA.lao2
#'  \item lrm - Passed to lm and lme
#'  \item llm - Passed to lm, lme and DA.llm
#'  \item llm2 - Passed to lm, lme and DA.llm2
#'  \item rai - Passed to raida
#'  \item spe - Passed to cor.test
#'  \item pea - Passed to cor.test
#' }
#' @return A list of results:
#' \itemize{
#'  \item table - Summary of results
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @export

allDA <- function(data, predictor, paired = NULL, tests = c("pea","spe","per","bay","adx","wil","ttt","ltt","ltt2","neb","erq","ere","msf","zig","ds2","lim","aov","lao","lao2","kru","lrm","llm","llm2","rai"), relative = TRUE, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list()){

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
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
  } else {
    count_table <- data
  }
  
  # Checks
  if(min(count_table) < 0 & is.numeric(predictor)) stop("Numeric predictor and negative values in count_table is currently not supported")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(levels(as.factor(predictor))) < 2) stop("predictor should have at least two levels")
  
  # Prune tests argument
  tests <- unique(tests)
  tests <- prune.tests.DA(tests, predictor, paired, relative)
  
  # Set seed
  set.seed(rng.seed)
  message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Run tests
  results <- run.tests.DA(count_table, predictor, paired, tests, relative, args, cores, p.adj)
  
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
  
  # Add taxa data
  if(class(data) == "phyloseq"){
    if(!is.null(tax_table(data, errorIfNULL = FALSE))){
      newresults <- list()
      tax <- tax_table(data)
      for(i in 1:length(results)){
        subres <- results[[i]] 
        subres <- merge(subres, tax, by.x = "Feature", by.y = "row.names")
        newresults[[i]] <- subres
      }
      names(newresults) <- names(results)
    
      df.combined <- merge(df.combined, tax, by.x = "Feature", by.y = "row.names")
      rownames(df.combined) <- NULL
      
    } else {
      newresults <- results
    } 
  } else {
    newresults <- results
  }
  
  
  return(list(table = df.combined,results = newresults))

}




