#' Comparing differential abundance methods by FPR and AUC
#'
#' Calculating false positive rates and AUC (Area Under the receiver operator Curve) for various differential abundance methods
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param covars Either a named list with covariates, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param R Integer. Number of times to run the tests. Default 10
#' @param tests Character. Which tests to include. Default all (See below for details)
#' @param relative Logical. Should abundances be made relative? Only for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea". Default TRUE
#' @param effectSize Integer. The effect size for the spike-ins. Default 5
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. k=c(5,10,15): 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default c(5,5,5)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.anova If TRUE (default) linear models will output results and p-values from anova/drop1. If FALSE will output results for 2. level of the predictor.
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
#'  \item erq - EdgeR - Quasi likelihood - TMM normalization
#'  \item ere - EdgeR - Exact test - TMM normalization
#'  \item erq2 - EdgeR - Quasi likelihood - RLE normalization
#'  \item ere2 - EdgeR - Exact test - RLE normalization
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
#'  \item pea - Pearson correlation
#'  \item poi - Poisson GLM with log of library size as offset
#'  \item qpo - Quasi-Poisson GLM with log of library size as offset
#'  \item vli - Limma with voom
#'  \item zpo - Zero-inflated Poisson GLM
#'  \item znb - Zero-inflated Negative Binomial GLM
#'  \item fri - Friedman Rank Sum test
#'  \item qua - Quade test
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
#'  \item erq - Passed to calcNormFactors, estimateDisp, glmQLFit and glmQLFTest
#'  \item ere - Passed to calcNormFactors, estimateCommonDisp, estimateTagwiseDisp and exactTest
#'  \item msf - Passed to fitFeatureModel
#'  \item zig - Passed to fitZig
#'  \item ds2 - Passed to DESeq
#'  \item lim - Passed to eBayes and lmFit
#'  \item lli - Passed to eBayes, lmFit and DA.lli
#'  \item lli2 - Passed to eBayes, lmFit and DA.lli2
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
#'  \item poi - Passed to glm/glmer
#'  \item qpo - Passed to glm
#'  \item vli - Passed to voom, eBayes and lmFit
#'  \item zpo - Passed to zeroinfl
#'  \item znb - Passed to zeroinfl
#'  \item fri - Passed to friedman.test
#'  \item qua - Passed to quade.test
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

testDA <- function(data, predictor, paired = NULL, covars = NULL, R = 10, tests = c("qua","fri","zpo","znb","vli","qpo","poi","pea","neb","rai","per","bay","adx","wil","ttt","ltt","ltt2","erq","erq2","ere","ere2","msf","zig","ds2","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe"), relative = TRUE, effectSize = 5, k = c(5,5,5), cores = (detectCores()-1), rng.seed = 123, args = list(), out.anova = TRUE){

  stopifnot(exists("data"))
  
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
    if(!is.null(covars)){
      covars.n <- covars
      covars <- list()
      for(i in 1:length(covars.n)){
        covars[[i]] <- suppressWarnings(as.matrix(sample_data(data)[,covars.n[i]]))
      }
      names(covars) <- covars.n
    } 
  } else {
    count_table <- data
  }
  
  # Coerce data
  if(!is.null(paired)) paired <- as.factor(paired)
  count_table <- as.matrix(count_table)
  
  # Checks
  if(min(count_table) < 0) stop("Count_table contains negative values!")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(levels(as.factor(predictor))) < 2) stop("predictor should have at least two levels")
  
  # Prune tests argument
  tests <- unique(tests)
  tests <- prune.tests.DA(tests, predictor, paired, covars, relative)
  tests.par <- paste0(unlist(lapply(1:R, function(x) rep(x,length(tests)))),"_",rep(tests,R))
  
  # neb warning
  if("neb" %in% tests & !is.null(paired)){
    message("As 'neb' is included and a 'paired' variable is supplied this might take a long time")
  }
  
  # Set seed
  set.seed(rng.seed)
  message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Numeric predictor
  if(is.numeric(predictor[1])){
    num.pred <- TRUE
    message("predictor is assumed to be a continuous/quantitative variable")
  } else {
    num.pred <- FALSE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      if(is.numeric(covars[[i]][1])){
        message(paste(names(covars)[i],"is assumed to be a continuous/quantitative variable"))
      } else {
        message(paste(names(covars)[i],"is assumed to be a categorical variable"))
      }
    }
  }
  
  # Shuffle predictor
  if(is.null(paired)){
    rands <- lapply(1:R,function(x) sample(predictor))
  } else {
    rands <- lapply(1:R,function(x) unsplit(lapply(split(predictor,paired), sample), paired))
  }
  
  # Spikeins
  spikeds <- lapply(1:R,function(x) spikein(count_table, rands[[x]], effectSize,  k, num.pred, relative))
  count_tables <- lapply(1:R,function(x) spikeds[[x]][[1]])
  
  ### Run tests
  # Progress bar
  pb <- txtProgressBar(max = length(tests.par), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Start parallel
  if(cores == 1) {
    registerDoSEQ() 
  } else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  }

  # Run the tests in parallel
  results <- foreach(i = tests.par , .options.snow = opts) %dopar% {

    # Extract run info
    run.no <- as.numeric(gsub("_.*","",i))
    i <- gsub(".*_","",i)
    
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
    
    # Run tests
    res.sub <- tryCatch(switch(i,
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative),wil.args)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative),ttt.args)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative),ltt.args)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired),ltt2.args)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.anova),neb.args)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars),erq.args)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),ere.args)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars),erq2.args)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),ere2.args)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),msf.args)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars),zig.args)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars),ds2.args)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative),per.args)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired),bay.args)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),adx.args)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.anova),lim.args)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.anova),lli.args)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.anova),lli2.args)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], relative),kru.args)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, relative),aov.args)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative),lao.args)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars),lao2.args)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars, relative,out.anova),lrm.args)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.anova),llm.args)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.anova),llm2.args)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),rai.args)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative),spe.args)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative),pea.args)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.anova),poi.args)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,out.anova),qpo.args)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.anova),vli.args)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,out.anova),zpo.args)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,out.anova),znb.args)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative),fri.args)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative),qua.args))),
                        
                        error = function(e) NULL)
    
    if(!is.null(res.sub)){
      res.sub[is.na(res.sub$pval),"pval"] <- 1
    }

    return(res.sub)
    
  }
  names(results) <- tests.par
  
  # Handle failed tests
  results <- results[!sapply(results,is.null)]
  
  if(length(unique(gsub(".*_","",names(results)))) != length(tests)){
    if(length(tests) - length(unique(gsub(".*_","",names(results)))) == 1){
      message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"was excluded due to failure"))
    } else {
      message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"were excluded due to failure"))
    }
    tests <- unique(gsub(".*_","",names(results)))
  }
  
  final.results <- foreach(r = 1:R) %do% {

    res.sub <- results[names(results)[gsub("_.*","",names(results)) == r]]
    
    # Split ALDEx2 results in t.test and wilcoxon
    if("adx" %in% tests){
      adx.t <- as.data.frame(res.sub[paste0(r,"_","adx")])[,c(1:7,12)]
      adx.w <- as.data.frame(res.sub[paste0(r,"_","adx")])[,c(1:7,12)]
      colnames(adx.t) <- gsub(".*_adx.","",colnames(adx.t))
      colnames(adx.w) <- colnames(adx.t)
      adx.t$pval <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,8])
      adx.w$pval <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,10])
      adx.t$pval.adj <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,12])
      adx.w$pval.adj <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,13])
      adx.t$Method <- "ALDEx2 t-test (adx)"
      adx.w$Method <- "ALDEx2 wilcox (adx)"
      res.sub[paste0(r,"_","adx")] <- NULL
      res.names <- names(res.sub)
      res.sub <- c(res.sub,list(adx.t),list(adx.w))
      names(res.sub) <- c(res.names,paste0(r,"_","adx.t"),paste0(r,"_","adx.w"))
    }

    # Insert spiked column
    newnames <- gsub(".*_","",names(res.sub))
    res.sub <- foreach(rsp = 1:length(res.sub)) %do% {
      temp <- res.sub[[rsp]]
      temp$Spiked <- "No"
      temp[temp$Feature %in% spikeds[[r]][[2]],"Spiked"] <- "Yes"
      return(temp)
    }
    names(res.sub) <- newnames
    
    # Confusion matrix
    totalPos <- sapply(res.sub,function(x) nrow(x[x$pval < 0.05,]))
    totalNeg <- sapply(res.sub,function(x) nrow(x[x$pval >= 0.05,])) 
    trueNeg <- totalNeg  #if effectSize == 1
    truePos <- 0  #if effectSize == 1
    falseNeg <- 0 #if effectSize == 1
    if(effectSize != 1){
      truePos <- sapply(res.sub, function(x) sum(x[x$pval < 0.05,"Feature"] %in% spikeds[[r]][[2]]))
      falseNeg <- sapply(res.sub, function(x) sum(x[x$pval >= 0.05,"Feature"] %in% spikeds[[r]][[2]]))
    }
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR 
    fprs <- sapply(1:length(res.sub), function(x) {
      if(totalPos[x] != 0){
        falsePos[x] / (totalPos[x] + totalNeg[x])
      } else {0}})
    
    
    # Spike detection rate
    sdrs <- sapply(1:length(res.sub), function(x) truePos[x] / sum(k))
    
    # AUC
    aucs <- sapply(1:length(res.sub), function(x) {
      if(effectSize != 1){
        test_roc <- NULL
        tryCatch(
          test_roc <- pROC::roc(as.numeric(res.sub[[x]]$Feature %in% spikeds[[r]][[2]]) ~ res.sub[[x]]$pval, auc=TRUE, direction = ">"),
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
    df.combined <- data.frame(Method = sapply(res.sub, function(x) x$Method[1]),
                              AUC = aucs,
                              FPR = fprs,
                              Spike.detect.rate = sdrs,
                              Run = r)
    rownames(df.combined) <- NULL

    return(list(df.combined,res.sub))
    
  }
  
  output.results <- do.call(rbind,lapply(final.results, function(x) x[[1]]))
  output.all.results <- lapply(final.results, function(x) x[[2]])
  
  out <- list(table = output.results, results = output.all.results)
  class(out) <- "DA"
  return(out)
}




