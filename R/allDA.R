#' Run many differential abundance methods
#'
#' Run many differential abundance tests at a time
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2", "zig" and "fri"
#' @param covars Either a named list with covariates, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param tests Character. Which tests to include. Default all (Except ANCOM, see below for details)
#' @param relative Logical. If TRUE (default) abundances are made relative for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea", and there is an offset of log(LibrarySize) for "neb", "poi", "qpo", "zpo" and "znb"
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for pvalue adjustment. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.anova If TRUE linear models will output results and p-values from anova/drop1. If FALSE will output results for 2. level of the predictor. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param alpha P-value threshold for calling significance. Default 0.05
#' @param core.check If TRUE will make an interactive check that the amount of cores specified are desired. Only if cores>10. This is to ensure that the function doesn't automatically overloads a server with workers.  
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
#'  \item raw - A data.frame with raw p-values from all methods
#'  \item adj - A data.frame with adjusted p-values from all methods (detection/no-detection from anc and sam)
#'  \item est - A data.frame with estimates/fold.changes from all relevant methods
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: $results[[2]]["wil"]
#' }
#' 
#' @export

allDA <- function(data, predictor, paired = NULL, covars = NULL, tests = c("neb","per","bay","adx","sam","qua","fri","znb","zpo","vli","qpo","poi","pea","spe","wil","ttt","ltt","ltt2","erq","ere","erq2","ere2","msf","zig","ds2","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","rai"), relative = TRUE, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list(), out.anova = NULL, alpha = 0.05, core.check = TRUE){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(core.check){
    if(cores > 10){
      ANSWER <- readline(paste("You are about to run allDA using",cores,"cores. Enter y to proceed "))
      if(ANSWER != "y") stop("Process aborted")
    }
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
    predictor <- unlist(sample_data(data)[,predictor])
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    if(!is.null(covars)){
      covars.n <- covars
      covars <- list()
      for(i in 1:length(covars.n)){
        covars[[i]] <- unlist(sample_data(data)[,covars.n[i]])
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
  if(length(unique(predictor)) < 2) stop("predictor should have at least two levels")
  
  # Prune tests argument
  tests <- unique(tests)
  if(!"zzz" %in% tests) tests <- prune.tests.DA(tests, predictor, paired, covars, relative)
  if(length(tests) == 0) stop("No tests to run!")
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # predictor
  if(is.numeric(predictor)){
    message("predictor is assumed to be a quantitative variable")
    if(length(levels(as.factor(predictor))) == 2){
      ANSWER <- readline("The predictor is quantitative, but only contains 2 unique values. Are you sure this is correct? Enter y to proceed ")
      if(ANSWER != "y") stop("Wrap the predictor with as.factor(predictor) to treat it is a categorical variable")
    }
  } else {
    if(length(levels(as.factor(predictor))) > length(unique(predictor))) stop("predictor has more levels than unique values!")
    message(paste("predictor is assumed to be a categorical variable with",length(unique(predictor)),"levels:",paste(levels(as.factor(predictor)),collapse = ", ")))
  }

  # Out.anova
  if(is.null(out.anova)){
    if(is.numeric(predictor)) out.anova <- FALSE
    if(length(unique(predictor)) == 2) out.anova <- FALSE
    if(length(unique(predictor)) > 2) out.anova <- TRUE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      if(is.numeric(covars[[i]][1])){
        message(paste(names(covars)[i],"is assumed to be a quantitative variable"))
      } else {
        message(paste(names(covars)[i],"is assumed to be a categorical variable with",length(unique(covars[[i]])),"levels:",paste(levels(as.factor(covars[[i]])),collapse = ", ")))
      }
    }
  }
  
  # Set seed
  message(paste("Seed is set to",rng.seed))
  set.seed(rng.seed)
  
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

    # Set seed
    set.seed(rng.seed)
    
    # Extract test arguments
    if(!all(names(args) %in% tests)) stop("One or more names in list with additional arguments does not match names of tests")
    for(j in seq_along(args)){
      assign(paste0(names(args)[j],".DAargs"),args[[j]],pos=1)
    }
    test.args <- paste0(tests,".DAargs")
    test.boo <- lapply(test.args,exists)
    for(l in seq_along(test.args)){
      if(test.boo[l] == FALSE) assign(test.args[l], list(),pos=1)
    }
    
    if(!is.na(pmatch("zzz",i))){
      zzz.DAargs <- get(paste0(i,".DAargs"))
      i <- "zzz"
    } 
    
    on.exit(suppressWarnings(rm(list=test.args, pos = 1)), add = TRUE)
    
    res.sub <- tryCatch(switch(i,
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),zzz.DAargs)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),wil.DAargs)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),ttt.DAargs)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative, p.adj),ltt.DAargs)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, p.adj),ltt2.DAargs)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),neb.DAargs)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),erq.DAargs)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere.DAargs)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),erq2.DAargs)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere2.DAargs)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),msf.DAargs)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),zig.DAargs)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),ds2.DAargs)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),per.DAargs)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),bay.DAargs)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),adx.DAargs)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),lim.DAargs)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),lli.DAargs)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),lli2.DAargs)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, relative, p.adj),kru.DAargs)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, relative, p.adj),aov.DAargs)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative, p.adj),lao.DAargs)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, p.adj),lao2.DAargs)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, relative,out.anova, p.adj),lrm.DAargs)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),llm.DAargs)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),llm2.DAargs)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),rai.DAargs)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),spe.DAargs)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),pea.DAargs)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.anova, p.adj),poi.DAargs)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),qpo.DAargs)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.anova, p.adj),vli.DAargs)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),zpo.DAargs)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.anova, p.adj),znb.DAargs)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),fri.DAargs)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),qua.DAargs)),
                               anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),anc.DAargs)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),sam.DAargs))),
                        
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
  
  # Raw p-values
  Pval.raw <- lapply(results,function(x) tryCatch(as.data.frame(x[,c("Feature","pval")]), error = function(e) NULL))
  Pval.raw <- Pval.raw[!sapply(Pval.raw,is.null)]
  if(length(Pval.raw) > 0){
    df.raw <- suppressWarnings(Reduce(function(x,y) merge(x, y, by= "Feature", all.x = TRUE, all.y = TRUE), Pval.raw))
    colnames(df.raw)[2:ncol(df.raw)] <- names(Pval.raw)
    if(class(data) == "phyloseq") df.raw <- add.tax.DA(data, df.raw)
  } else {
    df.raw <- NULL
  }

  # Adjusted p-values
  Pval.adj <- lapply(results,function(x) tryCatch(as.data.frame(x[,c("Feature","pval.adj")]), error = function(e) NULL))
  Pval.adj <- Pval.adj[!sapply(Pval.adj,is.null)]
  if(length(Pval.adj) > 0){
    df.adj <- suppressWarnings(Reduce(function(x,y) merge(x, y, by= "Feature", all.x = TRUE, all.y = TRUE), Pval.adj))
    colnames(df.adj)[2:ncol(df.adj)] <- names(Pval.adj)
    if("sam" %in% names(results)) {
      if("Sig" %in% colnames(results$sam)){
        df.adj <- merge(df.adj, results$sam[,c("Feature","Sig")], by = "Feature")
      } else {
        sam.adj <- results$sam[,c("Feature","Sig.up","Sig.lo")]
        sam.adj$Sig <- "No"
        sam.adj[sam.adj$Sig.up == "Yes","Sig"] <- "Up"
        sam.adj[sam.adj$Sig.lo == "Yes","Sig"] <- "Down"
        df.adj <- merge(df.adj, sam.adj[,c("Feature","Sig")], by = "Feature")  
      }
      colnames(df.adj)[ncol(df.adj)] <- "sam"
    }
    if("anc" %in% names(results)) {
      df.adj <- merge(df.adj, results$anc[,c("Feature","Detected")], by = "Feature")
      colnames(df.adj)[ncol(df.adj)] <- "anc"
    }
    if(class(data) == "phyloseq") df.adj <- add.tax.DA(data, df.adj)
  } else {
    df.adj <- NULL
  }
  
  ## Estimate
  est.name <- list(sam = "Fold.change",
                   znb = "Estimate",
                   zpo = "Estimate",
                   qpo = "Estimate",
                   poi = "Estimate",
                   neb = "Estimate",
                   lrm = "Estimate",
                   llm = "Estimate",
                   llm2 = "Estimate",
                   vli = c("logFC",paste0("predictor",levels(as.factor(predictor))[2])),
                   lim = c("logFC",paste0("predictor",levels(as.factor(predictor))[2])),
                   lli = c("logFC",paste0("predictor",levels(as.factor(predictor))[2])),
                   lli2 = c("logFC",paste0("predictor",levels(as.factor(predictor))[2])),
                   pea = "cor",
                   spe = "rho",
                   per = "FC",
                   bay = "ordering",
                   adx.t = "effect",
                   adx.w = "effect",
                   wil = "FC",
                   ttt = "FC",
                   ltt = "FC",
                   ltt2 = "FC",
                   erq = c("logFC",paste0("logFC.predictor",levels(as.factor(predictor))[2])),
                   ere = "logFC",
                   erq2 = c("logFC",paste0("logFC.predictor",levels(as.factor(predictor))[2])),
                   ere2 = "logFC",
                   msf = "logFC",
                   zig = paste0("predictor",levels(as.factor(predictor))[2]),
                   ds2 = "log2FoldChange")

  if(!is.numeric(predictor) & length(unique(predictor)) > 2){
    df.est <- NULL
  } else {
    list.est <- foreach(ll = names(results)) %do% {
      if(ll %in% names(est.name)){
        if(any(est.name[ll][[1]] %in% colnames(results[ll][[1]]))){
          est.sub.name <- est.name[ll][[1]][est.name[ll][[1]] %in% colnames(results[ll][[1]])]
          est.sub <- results[ll][[1]][,c("Feature",est.sub.name)]
          colnames(est.sub) <- c("Feature",paste0(ll,"_",est.sub.name,""))
          return(est.sub)
        } else return(NULL)
      } else return(NULL)
    }
    list.est <- list.est[!sapply(list.est, is.null)]
    if(length(list.est) > 0){
      df.est <- Reduce(function(x,y) merge(x, y, by= "Feature", all.x = TRUE, all.y = TRUE), list.est)
      if(class(data) == "phyloseq") df.est <- add.tax.DA(data, df.est)
    } else {
      df.est <- NULL
    } 
  }

  # Add tax table to results
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
    } else {
      newresults <- results
    } 
  } else {
    newresults <- results
  }
  
  return(list(raw = df.raw, adj = df.adj, est = df.est,results = newresults))

}




