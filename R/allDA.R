#' Run many differential abundance methods
#'
#' Run many differential abundance tests at a time
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2", "zig" and "fri"
#' @param covars Either a named list with covariates, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param tests Character. Which tests to include. Default all (Except ANCOM, see below for details)
#' @param relative Logical. Should abundances be made relative? Only for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea". Default TRUE
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for pvalue adjustment. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.anova If TRUE (default) linear models will output results and p-values from anova/drop1. If FALSE will output results for 2. level of the predictor.
#' @param alpha P-value threshold for calling significance. Default 0.05
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
#'  \item anc - ANCOM (by default not included, as it is very slow)
#'  \item sam - SAMSeq
#'  \item zzz - A user-defined method (See ?DA.zzz)
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
#'  \item anc - Passed to ANCOM
#'  \item sam - Passed to SAMseq
#' }
#' @return A list of results:
#' \itemize{
#'  \item table - Summary of results
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @export

allDA <- function(data, predictor, paired = NULL, covars = NULL, tests = c("sam","qua","fri","znb","zpo","vli","qpo","poi","pea","spe","per","bay","adx","wil","ttt","ltt","ltt2","neb","erq","ere","erq2","ere2","msf","zig","ds2","lim","aov","lao","lao2","kru","lrm","llm","llm2","rai"), relative = TRUE, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list(), out.anova = TRUE, alpha = 0.05){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(cores > 10){
    ANSWER <- readline(paste("You are about to run allDA using",cores,"cores. Enter y to proceed "))
    if(ANSWER != "y") stop("Process aborted")
  }
  
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

  # Checks
  if(relative) if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) stop("Count_table must only contain integer values")
  if(min(count_table) < 0) stop("Count_table contains negative values!")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(levels(as.factor(predictor))) < 2) stop("predictor should have at least two levels")
  
  # Prune tests argument
  tests <- unique(tests)
  if(!"zzz" %in% tests) tests <- prune.tests.DA(tests, predictor, paired, covars, relative)
  
  # Set seed
  set.seed(rng.seed)
  message(paste("Seed is set to",rng.seed))
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Numeric predictor
  if(is.numeric(predictor[1])){
    message("predictor is assumed to be a continuous/quantitative variable")
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
  
  # Run tests
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
    on.exit(stopCluster(cl))
  }
  
  # Run tests in parallel
  results <- foreach(i = tests, .options.snow = opts) %dopar% {
    
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
    
    if(!is.na(pmatch("zzz",i))){
      zzz.args <- get(paste0(i,".args"))
      i <- "zzz"
    } 
    
    res.sub <- tryCatch(switch(i,
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),zzz.args)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),wil.args)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),ttt.args)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative, p.adj),ltt.args)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, p.adj),ltt2.args)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),neb.args)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),erq.args)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere.args)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),erq2.args)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere2.args)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),msf.args)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),zig.args)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),ds2.args)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),per.args)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),bay.args)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),adx.args)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),lim.args)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),lli.args)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),lli2.args)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, relative, p.adj),kru.args)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, relative, p.adj),aov.args)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative, p.adj),lao.args)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, p.adj),lao2.args)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, relative,out.anova, p.adj),lrm.args)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),llm.args)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),llm2.args)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),rai.args)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),spe.args)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),pea.args)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),poi.args)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),qpo.args)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),vli.args)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),zpo.args)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),znb.args)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),fri.args)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),qua.args)),
                               anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),anc.args)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),sam.args))),
                        
                        error = function(e) NULL)
    
    if(!is.null(res.sub) & !i %in% c("sam","anc")){
      res.sub[is.na(res.sub$pval),"pval"] <- 1
      res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
    }
    
    return(res.sub)
    
  }
  names(results) <- tests
  
  # Handle failed tests
  results <- results[!sapply(results,is.null)]
  if(length(names(results)) != length(tests)){
    if(length(tests) - length(names(results)) == 1){
      message(paste(paste(tests[!tests %in% names(results)],collapse = ", "),"was excluded due to failure"))
    } else {
      message(paste(paste(tests[!tests %in% names(results)],collapse = ", "),"were excluded due to failure"))
    }
    tests <- names(results)
  }
  
  # Split ALDEx2 results in t.test and wilcoxon
  if("adx" %in% names(results)){
    adx.t <- as.data.frame(results["adx"])[,c(1:7,14)]
    adx.w <- as.data.frame(results["adx"])[,c(1:7,14)]
    colnames(adx.t) <- gsub("adx.","",colnames(adx.t))
    colnames(adx.w) <- colnames(adx.t)
    adx.t$pval <- as.numeric(as.data.frame(results["adx"])$adx.we.ep)
    adx.w$pval <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep)
    adx.t$pval.adj <- as.numeric(as.data.frame(results["adx"])$adx.we.ep.adj)
    adx.w$pval.adj <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep.adj)
    adx.t$Method <- "ALDEx2 t-test (adx)"
    adx.w$Method <- "ALDEx2 wilcox (adx)"
    results["adx"] <- NULL
    res.names <- names(results)
    results <- c(results,list(adx.t),list(adx.w))
    names(results) <- c(res.names,"adx.t","adx.w")
  }

  # P-values for ANCOM and SAMseq
  if("anc" %in% names(results)){
    ancdf <- as.data.frame(results[["anc"]])
    ancdf$pval.adj <- 1
    ancdf[ancdf$Detected == "Yes","pval.adj"] <- 0
    ancdf$pval <- NA
    results["anc"] <- NULL
    res.names <- names(results)
    results <- c(results,list(ancdf))
    names(results) <- c(res.names,"anc")
  }
  if("sam" %in% names(results)){
    samdf <- as.data.frame(results[["sam"]])
    samdf$pval.adj <- 1
    if("Sig" %in% colnames(samdf)){
      samdf[samdf$Sig == "Yes","pval.adj"] <- 0
    }
    if("Sig.up" %in% colnames(samdf)){
      samdf[samdf$Sig.up == "Yes","pval.adj"] <- 0
    }
    if("Sig.lo" %in% colnames(samdf)){
      samdf[samdf$Sig.lo == "Yes","pval.adj"] <- 0
    }
    samdf$pval <- NA
    results["sam"] <- NULL
    res.names <- names(results)
    results <- c(results,list(samdf))
    names(results) <- c(res.names,"sam")
  }
  
  # Positives
  Pos.raw <- sapply(results,function(x) x[x$pval < alpha,"Feature"])
  Pos.adj <- sapply(results,function(x) x[x$pval.adj < alpha,"Feature"])

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




