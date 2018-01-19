#' Comparing differential abundance/expression methods by FPR and AUC
#'
#' Calculating false positive rates and AUC (Area Under the Receiver Operating Characteristic (ROC) Curve) for various differential abundance and expression methods
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. Only for "anc", "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "ds2x", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param R Integer. Number of times to run the tests. Default 10
#' @param tests Character. Which tests to include. Default all (Except ANCOM and mvabund, see below for details)
#' @param relative Logical. If TRUE (default) abundances are made relative for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea", and there is an offset of \code{log(LibrarySize)} for "neb", "poi", "qpo", "zpo" and "znb"
#' @param effectSize Numeric. The effect size for the spike-ins. Default 5
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. \code{k=c(5,10,15)}: 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default NULL, which will spike 2 percent of the total amount of features in each tertile (a total of 6 percent), but minimum c(5,5,5)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.all If TRUE linear models will output results and p-values from \code{anova}/\code{drop1}, ds2/ds2x will run LRT and not Wald test, erq and erq2 will produce one p-value for the predictor, and lim, lli, lli2, lim, vli will run F-tests. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param alpha q-value threshold for determining significance for \code{Spike.detect.rate}. Default 0.1
#' @param core.check If TRUE will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.  
#' @param verbose If TRUE will print informative messages
#' @details Currently implemented methods:
#' \itemize{
#'  \item per - Permutation test with user defined test statistic
#'  \item bay - baySeq
#'  \item adx - ALDEx t-test and wilcoxon
#'  \item wil - Wilcoxon Rank Sum on relative abundances
#'  \item ttt - Welch t.test on relative abundances
#'  \item ltt - Welch t.test, but reads are first transformed with \code{log(abundance + delta1)} then turned into relative abundances
#'  \item ltt2 - Welch t.test, but with relative abundances transformed with \code{log(relative abundance + delta2)}
#'  \item neb - Negative binomial GLM with log of library size as offset
#'  \item erq - EdgeR - Quasi likelihood - TMM normalization
#'  \item ere - EdgeR - Exact test - TMM normalization
#'  \item erq2 - EdgeR - Quasi likelihood - RLE normalization
#'  \item ere2 - EdgeR - Exact test - RLE normalization
#'  \item msf - MetagenomeSeq feature model
#'  \item zig - MetagenomeSeq zero-inflated gaussian
#'  \item ds2 - DESeq2
#'  \item ds2x - DESeq2 with manual geometric means
#'  \item lim - LIMMA. Moderated linear models based on emperical bayes
#'  \item lli - LIMMA, but reads are first transformed with \code{log(abundance + delta1)} then turned into relative abundances
#'  \item lli2 - LIMMA, but with relative abundances transformed with \code{log(relative abundance + delta2)}
#'  \item kru - Kruskal-Wallis on relative abundances
#'  \item aov - ANOVA on relative abundances
#'  \item lao - ANOVA, but reads are first transformed with \code{log(abundance + delta1)} then turned into relative abundances
#'  \item lao2 - ANOVA, but with relative abundances transformed with \code{log(relative abundance + delta2)}
#'  \item lrm - Linear regression on relative abundances
#'  \item llm - Linear regression, but reads are first transformed with \code{log(abundance + delta1)} then turned into relative abundances
#'  \item llm2 - Linear regression, but with relative abundances transformed with \code{log(relative abundance + delta2)}
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
#'  \item mva - mvabund (by default not included, as it is very slow)
#'  \item zzz - A user-defined method (See \code{?DA.zzz})
#' }
#' "neb" can be slow if there is a paired argument.
#' 
#' "anc" and "mva" are slow compared to the other methods.
#' 
#' Additional arguments can be passed to the internal functions with the \code{args} argument. 
#' It should be structured as a list with elements named by the tests: 
#' E.g. passing to the \code{DA.per} function that it should only run 1000 iterations: \code{args = list(per=list(noOfIterations=1000))}.
#' Include that the log t.test should use a pseudocount of 0.1: \code{args = list(per=list(noOfIterations=1000), ltt=list(delta=0.1))}. 
#' Additional arguments are simply seperated by commas.
#' 
#' Below is an overview of which functions get the arguments that are passed to a specific test
#' \itemize{
#'  \item per - Passed to \code{DA.per}
#'  \item bay - Passed to \code{getPriors}, \code{getLikelihoods} and \code{DA.bay}
#'  \item adx - Passed to \code{aldex} and \code{DA.adx}
#'  \item wil - Passed to \code{wilcox.test} and \code{DA.wil}
#'  \item ttt - Passed to \code{t.test} and \code{DA.ttt}
#'  \item ltt - Passed to \code{t.test} and \code{DA.ltt}
#'  \item ltt2 - Passed to \code{t.test} and \code{DA.ltt2}
#'  \item neb - Passed to \code{glm.nb}, \code{glmer.nb} and \code{DA.neb}
#'  \item erq(2) - Passed to \code{calcNormFactors}, \code{estimateDisp}, \code{glmQLFit}, \code{glmQLFTest} and \code{DA.erq}
#'  \item ere(2) - Passed to \code{calcNormFactors}, \code{estimateCommonDisp}, \code{estimateTagwiseDisp}, \code{exactTest} and \code{DA.ere}
#'  \item msf - Passed to \code{fitFeatureModel} and \code{DA.msf}
#'  \item zig - Passed to \code{fitZig} and \code{DA.zig}
#'  \item ds2(x) - Passed to \code{DESeq} and \code{DA.ds2}
#'  \item lim - Passed to \code{eBayes}, \code{lmFit} and \code{DA.lim}
#'  \item lli - Passed to \code{eBayes}, \code{lmFit} and \code{DA.lli}
#'  \item lli2 - Passed to \code{eBayes}, \code{lmFit} and \code{DA.lli2}
#'  \item kru - Passed to \code{kruskal.test} and \code{DA.kru}
#'  \item aov - Passed to \code{aov} and \code{DA.aov}
#'  \item lao - Passed to \code{aov} and \code{DA.lao}
#'  \item lao2 - Passed to \code{aov} and \code{DA.lao2}
#'  \item lrm - Passed to \code{lm}, \code{lme} and \code{DA.lrm}
#'  \item llm - Passed to \code{lm}, \code{lme} and \code{DA.llm}
#'  \item llm2 - Passed to \code{lm}, \code{lme} and \code{DA.llm2}
#'  \item rai - Passed to \code{raida} and \code{DA.rai}
#'  \item spe - Passed to \code{cor.test} and \code{DA.spe}
#'  \item pea - Passed to \code{cor.test} and \code{DA.pea}
#'  \item poi - Passed to \code{glm}, \code{glmer} and \code{DA.poi}
#'  \item qpo - Passed to \code{glm} and \code{DA.qpo}
#'  \item vli - Passed to \code{voom}, \code{eBayes}, \code{lmFit} and \code{DA.vli}
#'  \item zpo - Passed to \code{zeroinfl} and \code{DA.zpo}
#'  \item znb - Passed to \code{zeroinfl} and \code{DA.znb}
#'  \item fri - Passed to \code{friedman.test} and \code{DA.fri}
#'  \item qua - Passed to \code{quade.test} and \code{DA.qua}
#'  \item anc - Passed to \code{ANCOM} and \code{DA.anc}
#'  \item sam - Passed to \code{SAMseq} and \code{DA.sam}
#'  \item mva - Passed to \code{manyglm} and \code{summary.manyglm}
#' }
#' @return An object of class \code{DA}, which contains a list of results:
#' \itemize{
#'  \item table - FPR, AUC and spike detection rate for each run
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: \code{$results[[2]]["wil"]}
#'  \item details - A dataframe with details from the run
#'  \item run.times - A dataframe with average run times of the different methods
#' }
#' 
#' @import snow doSNOW foreach utils
#' @importFrom parallel detectCores
#' @importFrom pROC roc
#' @export

testDA <- function(data, predictor, paired = NULL, covars = NULL, R = 10, tests = c("neb","rai","per","bay","adx","sam","qua","fri","zpo","znb","vli","qpo","poi","pea","wil","ttt","ltt","ltt2","erq","erq2","ere","ere2","msf","zig","ds2","ds2x","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe"), relative = TRUE, effectSize = 5, k = NULL, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list(), out.all = NULL, alpha = 0.1, core.check = TRUE, verbose = TRUE){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(core.check){
    if(cores > 20){
      ANSWER <- readline(paste("You are about to run testDA using",cores,"cores. Enter y to proceed "))
      if(ANSWER != "y") stop("Process aborted")
    }
  }
  # Time taking
  t1 <- proc.time()
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    DAdata <- DA.phyloseq(data, predictor, paired, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
    covars <- DAdata$covars
  } else {
    count_table <- data
  }
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }
  
  # Coerce data
  if(!is.null(paired)) paired <- as.factor(paired)
  count_table <- as.matrix(count_table)
  
  # Checks
  if(relative) if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) stop("count_table must only contain integer values when relative=TRUE")
  if(min(count_table) < 0) stop("count_table contains negative values")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(unique(predictor)) < 2) stop("predictor should have at least two levels")

  # Prune tests argument
  decimal <- FALSE
  if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) decimal <- TRUE
  tests <- unique(tests)
  if(!"zzz" %in% tests) tests <- prune.tests.DA(tests, predictor, paired, covars, relative, decimal)
  tests.par <- paste0(unlist(lapply(1:R, function(x) rep(x,length(tests)))),"_",rep(tests,R))
  if(length(tests) == 0) stop("No tests to run!")
  
  # Run time warnings
  if(verbose){
    if("neb" %in% tests & !is.null(paired)){
      message("As 'neb' is included and a 'paired' variable is supplied, this might take a long time")
    } else {
      if("anc" %in% tests){
        message("As 'anc' is included, this might take some time")
      }
    }
  }

  # Set seed
  set.seed(rng.seed)
  if(verbose) message(paste("Seed is set to",rng.seed))

  # Create some random seeds for each run
  seeds <- rpois(R, lambda = (1:R)*1e6)

  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  if(nrow(count_table) <= 15) message("Warning: Dataset contains very few features") 
  
  # Spike vs no features
  if(is.null(k)){
    k <- rep(round(nrow(count_table)*0.02),3)
    if(sum(k) < 15){
      k <- c(5,5,5)
    } 
  } 
  if(sum(k) == nrow(count_table)) stop("Set to spike all features. Can't calculate FPR or AUC. Change k argument")
  if(sum(k) > nrow(count_table)) stop("Set to spike more features than are present in the data. Change k argument")
  if(sum(k) < 15 & sum(k) >= 10 & R <= 10) message("Few features spiked. Increase 'k' or set 'R' to more than 10 to ensure proper estimation of AUC and FPR")
  if(sum(k) < 10 & sum(k) >= 5 & R <= 20) message("Few features spiked. Increase 'k' or set 'R' to more than 20 to ensure proper estimation of AUC and FPR")                                  
  if(sum(k) < 5 & R <= 50) message("Very few features spiked. Increase 'k' or set 'R' to more than 50 to ensure proper estimation of AUC and FPR")
  if(sum(k) > nrow(count_table)/2) message("Set to spike more than half of the dataset, which might give unreliable estimates, Change k argument")   
                                 
  # predictor
  if(verbose) if(any(is.na(predictor))) message("Warning: Predictor contains NAs!")
  if(is.numeric(predictor[1])){
    num.pred <- TRUE
    if(verbose) message(paste("predictor is assumed to be a quantitative variable, ranging from",min(predictor, na.rm = TRUE),"to",max(predictor, na.rm = TRUE)))
    if(length(levels(as.factor(predictor))) == 2){
      ANSWER <- readline("The predictor is quantitative, but only contains 2 unique values. Are you sure this is correct? Enter y to proceed ")
      if(ANSWER != "y") stop("Wrap the predictor with as.factor(predictor) to treat it is a categorical variable")
    }
  } else {
    num.pred <- FALSE
    if(length(levels(as.factor(predictor))) > length(unique(predictor))) stop("predictor has more levels than unique values!")
    if(verbose) message(paste("predictor is assumed to be a categorical variable with",length(unique(predictor)),"levels:",paste(levels(as.factor(predictor)),collapse = ", ")))
  }
  
  # out.all
  if(is.null(out.all)){
    if(length(unique(predictor)) == 2) out.all <- FALSE
    if(length(unique(predictor)) > 2) out.all <- TRUE
    if(num.pred) out.all <- FALSE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      if(verbose) if(any(is.na(covars[[i]]))) message("Warning:",names(covars)[i],"contains NAs!")
      if(is.numeric(covars[[i]][1])){
        if(verbose) message(paste(names(covars)[i],"is assumed to be a quantitative variable, ranging from",min(covars[[i]], na.rm = TRUE),"to",max(covars[[i]], na.rm = TRUE)))
      } else {
        if(verbose) message(paste(names(covars)[i],"is assumed to be a categorical variable with",length(unique(covars[[i]])),"levels:",paste(levels(as.factor(covars[[i]])),collapse = ", ")))
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

    t1.sub <- proc.time()
    
    # Extract run info
    run.no <- as.numeric(gsub("_.*","",i))
    i <- gsub(".*_","",i)

    # Set subseed
    set.seed(seeds[run.no])
    
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
    
    # Run tests
    res.sub <- tryCatch(switch(i,
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars),zzz.DAargs)),
                               mva = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative, p.adj),mva.DAargs)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj),wil.DAargs)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj),ttt.DAargs)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj),ltt.DAargs)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, p.adj),ltt2.DAargs)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj),neb.DAargs)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),erq.DAargs)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj),ere.DAargs)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),erq2.DAargs)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj),ere2.DAargs)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj),msf.DAargs)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars, p.adj),zig.DAargs)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),ds2.DAargs)),
                               ds2x = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),ds2x.DAargs)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj),per.DAargs)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, p.adj),bay.DAargs)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]),adx.DAargs)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj),lim.DAargs)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj),lli.DAargs)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),lli2.DAargs)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], relative, p.adj),kru.DAargs)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, relative, p.adj),aov.DAargs)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative, p.adj),lao.DAargs)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, p.adj),lao2.DAargs)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars, relative,out.all, p.adj),lrm.DAargs)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj),llm.DAargs)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),llm2.DAargs)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj),rai.DAargs)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative, p.adj),spe.DAargs)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative, p.adj),pea.DAargs)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj),poi.DAargs)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj),qpo.DAargs)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj),vli.DAargs)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj),zpo.DAargs)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj),znb.DAargs)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj),fri.DAargs)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj),qua.DAargs)),
                               anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,sig = alpha),anc.DAargs)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,fdr.output = alpha),sam.DAargs))),
                        
                        error = function(e) NULL)
    
    if(!is.null(res.sub) & (!i %in% c("anc","sam","adx"))){
      res.sub[is.na(res.sub$pval),"pval"] <- 1
      res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
    }

    run.time.sub <- (proc.time()-t1.sub)[3]
    return(list(res.sub,run.time.sub))
    
  }
  run.times <- lapply(results, function(x) x[[2]])
  results <- lapply(results, function(x) x[[1]])
  
  names(results) <- tests.par
  names(run.times) <- tests.par

  # Handle failed tests
  results <- results[!sapply(results,is.null)]
  
  if(length(unique(gsub(".*_","",names(results)))) != length(tests)){
    if(length(tests) - length(unique(gsub(".*_","",names(results)))) == 1){
      message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"was excluded due to failure"))
    } else {
      message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"were excluded due to failure"))
    }

    # Produce informative messages
    if(all(tests[!tests %in% unique(gsub(".*_","",names(results)))] == "sam")){
      if(verbose) message("sam usually fails if some samples has too many zeroes")
    }
    if(all(c("sam","ere2","erq2","ds2x") %in% tests[!tests %in% unique(gsub(".*_","",names(results)))])){
      if(verbose) message("These tests usually fails if all features contain at least one zero")
    }
      
    tests <- unique(gsub(".*_","",names(results)))
  }
  
  final.results <- foreach(r = 1:R) %do% {

    res.sub <- results[names(results)[gsub("_.*","",names(results)) == r]]
    
    # Split ALDEx2 results in t.test and wilcoxon
    if("adx" %in% gsub(".*_","",names(res.sub))){
      adx.t <- as.data.frame(res.sub[paste0(r,"_","adx")])[,c(1:7,12:13)]
      adx.w <- as.data.frame(res.sub[paste0(r,"_","adx")])[,c(1:7,12:13)]
      colnames(adx.t) <- gsub(".*_adx.","",colnames(adx.t))
      colnames(adx.w) <- colnames(adx.t)
      adx.t$pval <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,8])
      adx.w$pval <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,10])
      adx.t$pval.adj <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,9])
      adx.w$pval.adj <- as.numeric(as.data.frame(res.sub[paste0(r,"_","adx")])[,11])
      adx.t$Method <- "ALDEx2 t-test (adx)"
      adx.w$Method <- "ALDEx2 wilcox (adx)"
      res.sub[paste0(r,"_","adx")] <- NULL
      res.names <- names(res.sub)
      res.sub <- c(res.sub,list(adx.t),list(adx.w))
      names(res.sub) <- c(res.names,paste0(r,"_","adx.t"),paste0(r,"_","adx.w"))
    }

    # Make pseudo-pval for ANCOM
    if("anc" %in% gsub(".*_","",names(res.sub))){
      anc <- as.data.frame(res.sub[paste0(r,"_","anc")])
      colnames(anc) <- gsub(".*_anc.","",colnames(anc))
      suppressWarnings(min.d <- min(anc[anc$Detected == "Yes","W"]))
      if(min.d == Inf) min.d <- max(anc$W)+1
      anc$pval <- 1/(anc$W+1) * 0.05/(1/(min.d+1))
      anc$pval.adj <- anc$pval
      res.sub[paste0(r,"_","anc")] <- NULL
      res.names <- names(res.sub)
      res.sub <- c(res.sub,list(anc))
      names(res.sub) <- c(res.names,paste0(r,"_","anc"))
    }
    
    # Make pseudo-pval for SAMseq
    if("sam" %in% gsub(".*_","",names(res.sub))){
      samdf <- as.data.frame(res.sub[paste0(r,"_","sam")])
      colnames(samdf) <- gsub(".*_sam.","",colnames(samdf))
      
      samdf$pval <- 1/rank(samdf$Score)
      samdf$pval.adj <- 1
      samdf[samdf$Sig == "Yes","pval.adj"] <- 0 
      
      res.sub[paste0(r,"_","sam")] <- NULL
      res.names <- names(res.sub)
      res.sub <- c(res.sub,list(samdf))
      names(res.sub) <- c(res.names,paste0(r,"_","sam"))
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
    totalPos <- sapply(res.sub,function(x) nrow(x[x$pval <= 0.05,]))
    totalNeg <- sapply(res.sub,function(x) nrow(x[x$pval > 0.05,])) 
    trueNeg <- totalNeg  #if effectSize == 1
    truePos <- 0  #if effectSize == 1
    falseNeg <- 0 #if effectSize == 1
    if(effectSize != 1){
      truePos <- sapply(res.sub, function(x) sum(x[x$pval <= 0.05,"Feature"] %in% spikeds[[r]][[2]]))
      falseNeg <- sapply(res.sub, function(x) sum(x[x$pval > 0.05,"Feature"] %in% spikeds[[r]][[2]]))
    }
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR 
    fprs <- sapply(1:length(res.sub), function(x) {
      if((falsePos[x] + trueNeg[x]) != 0){
        falsePos[x] / (falsePos[x] + trueNeg[x])
      } else {0}})
    
    
    # Spike detection rate
    # True positive for adjusted p-values
    truePos.adj <- 0  #if effectSize == 1
    if(effectSize != 1){
      truePos.adj <- sapply(res.sub, function(x) sum(x[x$pval.adj <= alpha,"Feature"] %in% spikeds[[r]][[2]]))
    }
    sdrs <- sapply(1:length(res.sub), function(x) truePos.adj[x] / sum(k))
    
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

    # SAMseq FDR
    if("sam" %in% newnames){
      if(nrow(res.sub[["sam"]][res.sub[["sam"]]$pval.adj <= alpha,]) == 0){
        df.combined[df.combined$Method == "SAMseq (sam)","FPR"] <- 0
      } else {
        df.combined[df.combined$Method == "SAMseq (sam)","FPR"] <- nrow(res.sub[["sam"]][res.sub[["sam"]]$pval.adj <= alpha & res.sub[["sam"]]$Spiked == "No",])/nrow(res.sub[["sam"]][res.sub[["sam"]]$pval.adj <= alpha,])
      }
    }
    
    return(list(df.combined,res.sub))
    
  }
  
  output.results <- do.call(rbind,lapply(final.results, function(x) x[[1]]))
  output.all.results <- lapply(final.results, function(x) x[[2]])
  
  if(num.pred){
    pred.det <- "Quantitative"
    pred.ord <- paste(min(predictor,na.rm=T),"to",max(predictor,na.rm=T))
  } else {
    if(length(levels(as.factor(predictor))) == 2){
      pred.det <- "Two-class"
      pred.ord <- paste(levels(as.factor(predictor))[1],"<",levels(as.factor(predictor))[2])
    } else {
      pred.det <- "Multi-class" 
      pred.ord <- paste(levels(as.factor(predictor)), collapse = ", ")
    }
  }
  
  if(is.null(paired)) pair.det <- "No" else pair.det <- "Yes"

  run.secs <- (proc.time() - t1)[3]

  if((run.secs)/60/60 > 1){
    run.time <- paste(round((run.secs)/60/60,2),"Hours")
   } else {
    run.time <- paste(round((run.secs)/60,2),"Minutes")
   }

  output.details <- data.frame(Features = nrow(count_table),
                               Samples = ncol(count_table),
                               Predictor = pred.det,
                               Ordering = pred.ord,
                               Paired = pair.det,
                               Covars = paste(names(covars), collapse = ", "),
                               RunTime = run.time,
                               Relative = relative,
                               EffectSize = effectSize,
                               Spiked = paste(paste0(c("Low:","Mid:","High:"),k), collapse = ", "),
                               RandomSeed = rng.seed,
                               OutAll = out.all)
  rownames(output.details) <- ""
  output.details <- as.data.frame(t(output.details))
  colnames(output.details) <- ""
  
  # Run times
  run.times.all <- foreach(i = unique(gsub(".*_","",names(results))),.combine = rbind) %do% {
    round(mean(as.numeric(run.times[gsub(".*_","",names(run.times)) == i]))/60,4)
  }
  run.times.all <- as.data.frame(run.times.all)
  rownames(run.times.all) <- unique(gsub(".*_","",names(results)))
  colnames(run.times.all) <- "Minutes"
  
  out <- list(table = output.results, results = output.all.results, details = output.details, run.times = run.times.all)
  class(out) <- "DA"
  
  # Create warning message
  output.summary.auc <- aggregate(AUC ~ Method, data = output.results, FUN = median)
  output.summary.fpr <- aggregate(FPR ~ Method, data = output.results, FUN = median)
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = output.results, FUN = median)
  df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
  df <- df[df$FPR <= 0.05,]
  if(nrow(df) > 0){
    df <- df[order(df$AUC, decreasing = TRUE),]
    if(df[1,"Spike.detect.rate"] == 0 & verbose) message("Spike detect rate is zero for the top method! You might run to re-run the analysis with a pruned dataset (see preDA) or a higher effectSize")
  }
  
  return(out)
}
