#' Run many differential abundance/expression methods
#'
#' Run many differential abundance and expression tests at a time, to easily compare their results
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "ds2x", "lrm", "llm", "llm2", "lim", "lli", "lli2", "zig" and "fri"
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param tests Character. Which tests to include. Default all (Except ANCOM, see below for details)
#' @param relative Logical. If TRUE (default) abundances are made relative for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea", and there is an offset of \code{log(LibrarySize)} for "neb", "poi", "qpo", "zpo" and "znb"
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param p.adj Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.all If TRUE models will output results and p-values from \code{anova}/\code{drop1}. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param alpha P-value threshold for calling significance. Default 0.05
#' @param core.check If TRUE (default) will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.  
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
#'  \item zzz - A user-defined method (See \code{?DA.zzz})
#' }
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
#' }
#' @return A list of results:
#' \itemize{
#'  \item raw - A data.frame with raw p-values from all methods
#'  \item adj - A data.frame with adjusted p-values from all methods (detection/no-detection from anc and sam)
#'  \item est - A data.frame with estimates/fold.changes from all relevant methods
#'  \item details - A dataframe with details from the run
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: \code{$results[[2]]["wil"]}
#' }
#' 
#' @export

allDA <- function(data, predictor, paired = NULL, covars = NULL, tests = c("neb","per","bay","adx","sam","qua","fri","znb","zpo","vli","qpo","poi","pea","spe","wil","ttt","ltt","ltt2","erq","ere","erq2","ere2","msf","zig","ds2","ds2x","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","rai"), relative = TRUE, cores = (detectCores()-1), rng.seed = 123, p.adj = "fdr", args = list(), out.all = NULL, alpha = 0.05, core.check = TRUE){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(core.check){
    if(cores > 20){
      ANSWER <- readline(paste("You are about to run allDA using",cores,"cores. Enter y to proceed "))
      if(ANSWER != "y") stop("Process aborted")
    }
  }
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    data <- DA.phyloseq(data, predictor, paired, covars)
    count_table <- data$count_table
    predictor <- data$predictor
    paired <- data$paired
    covars <- data$covars
  } else {
    count_table <- data
  }

  # Checks
  if(relative) if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) stop("Count_table must only contain integer values when relative=TRUE")
  if(min(count_table) < 0) stop("Count_table contains negative values!")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty!")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(unique(predictor)) < 2) stop("predictor should have at least two levels")
  
  # Prune tests argument
  decimal <- FALSE
  if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) decimal <- TRUE
  tests <- unique(tests)
  if(!"zzz" %in% tests) tests <- prune.tests.DA(tests, predictor, paired, covars, relative, decimal)
  if(length(tests) == 0) stop("No tests to run!")
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # predictor
  if(any(is.na(predictor))) message("Warning: Predictor contains NAs!")
  if(is.numeric(predictor)){
    message(paste("predictor is assumed to be a quantitative variable, ranging from",min(predictor, na.rm = TRUE),"to",max(predictor, na.rm = TRUE)))
    if(length(levels(as.factor(predictor))) == 2){
      ANSWER <- readline("The predictor is quantitative, but only contains 2 unique values. Are you sure this is correct? Enter y to proceed ")
      if(ANSWER != "y") stop("Wrap the predictor with as.factor(predictor) to treat it is a categorical variable")
    }
  } else {
    if(length(levels(as.factor(predictor))) > length(unique(predictor))) stop("predictor has more levels than unique values!")
    message(paste("predictor is assumed to be a categorical variable with",length(unique(predictor)),"levels:",paste(levels(as.factor(predictor)),collapse = ", ")))
  }

  # out.all
  if(is.null(out.all)){
    if(length(unique(predictor)) == 2) out.all <- FALSE
    if(length(unique(predictor)) > 2) out.all <- TRUE
    if(is.numeric(predictor)) out.all <- FALSE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in 1:length(covars)){
      if(any(is.na(covars[[i]]))) message("Warning:",names(covars)[i],"contains NAs!")
      if(is.numeric(covars[[i]][1])){
        message(paste(names(covars)[i],"is assumed to be a quantitative variable, ranging from",min(covars[[i]], na.rm = TRUE),"to",max(covars[[i]], na.rm = TRUE)))
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
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars),zzz.DAargs)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),wil.DAargs)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),ttt.DAargs)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative, p.adj),ltt.DAargs)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, p.adj),ltt2.DAargs)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.all, p.adj),neb.DAargs)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, out.all, p.adj),erq.DAargs)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere.DAargs)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, out.all, p.adj),erq2.DAargs)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),ere2.DAargs)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),msf.DAargs)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, p.adj),zig.DAargs)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.all, p.adj),ds2.DAargs)),
                               ds2x = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.all, p.adj),ds2x.DAargs)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired, relative, p.adj),per.DAargs)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),bay.DAargs)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor),adx.DAargs)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.all, p.adj),lim.DAargs)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.all, p.adj),lli.DAargs)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.all, p.adj),lli2.DAargs)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, relative, p.adj),kru.DAargs)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, relative, p.adj),aov.DAargs)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative, p.adj),lao.DAargs)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars, p.adj),lao2.DAargs)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars, relative,out.all, p.adj),lrm.DAargs)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.all, p.adj),llm.DAargs)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.all, p.adj),llm2.DAargs)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor, p.adj),rai.DAargs)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),spe.DAargs)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,relative, p.adj),pea.DAargs)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,relative,out.all, p.adj),poi.DAargs)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.all, p.adj),qpo.DAargs)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,covars,out.all, p.adj),vli.DAargs)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.all, p.adj),zpo.DAargs)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,covars,relative,out.all, p.adj),znb.DAargs)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),fri.DAargs)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired,relative,p.adj),qua.DAargs)),
                               anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),anc.DAargs)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,predictor,paired),sam.DAargs))),
                        
                        error = function(e) NULL)
    
    if(!is.null(res.sub) & !i %in% c("sam","anc","adx")){
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
    
    # Produce informative messages
    if(all(tests[!tests %in% unique(gsub(".*_","",names(results)))] == "sam")){
      message("sam usually fails if some samples has too many zeroes")
    }
    if(all(c("sam","ere2","erq2","ds2x") %in% tests[!tests %in% unique(gsub(".*_","",names(results)))])){
      message("sam, ere2, erq2 and ds2x usually fails if all features contain at least one zero")
    }
    
    tests <- names(results)
  }
  
  # Split ALDEx2 results in t.test and wilcoxon
  if("adx" %in% names(results)){
    adx.t <- as.data.frame(results["adx"])[,c(1:7,12:13)]
    adx.w <- as.data.frame(results["adx"])[,c(1:7,12:13)]
    colnames(adx.t) <- gsub("adx.","",colnames(adx.t))
    colnames(adx.w) <- colnames(adx.t)
    adx.t$pval <- as.numeric(as.data.frame(results["adx"])$adx.we.ep)
    adx.w$pval <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep)
    adx.t$pval.adj <- as.numeric(as.data.frame(results["adx"])$adx.we.eBH)
    adx.w$pval.adj <- as.numeric(as.data.frame(results["adx"])$adx.wi.eBH)
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
  est.name <- list(sam = "log2FC",
                   znb = "log2FC",
                   zpo = "log2FC",
                   qpo = "log2FC",
                   poi = "log2FC",
                   neb = "log2FC",
                   lrm = "log2FC",
                   llm = "log2FC",
                   llm2 = "log2FC",
                   vli = "logFC",
                   lim = "logFC",
                   lli = "logFC",
                   lli2 = "logFC",
                   pea = "cor",
                   spe = "rho",
                   per = "log2FC",
                   bay = "ordering",
                   adx.t = "effect",
                   adx.w = "effect",
                   wil = "log2FC",
                   ttt = "log2FC",
                   ltt = "log2FC",
                   ltt2 = "log2FC",
                   erq = c("logFC",paste0("logFC.predictor",levels(as.factor(predictor))[2])),
                   ere = "logFC",
                   erq2 = c("logFC",paste0("logFC.predictor",levels(as.factor(predictor))[2])),
                   ere2 = "logFC",
                   msf = "logFC",
                   zig = paste0("predictor",levels(as.factor(predictor))[2]),
                   ds2 = "log2FoldChange",
                   ds2x = "log2FoldChange",
                   rai = "log2FC")

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
  
  # Details
  if(is.numeric(predictor)){
    pred.det <- "Quantitative"
    pred.ord <- paste(min(predictor,na.rm=T),"to",max(predictor,na.rm=T))
  } else {
    if(length(levels(as.factor(predictor))) == 2){
      pred.det <- "Two-class"
      pred.ord <- paste0(levels(as.factor(predictor))[2],">",levels(as.factor(predictor))[1])
    } else {
      pred.det <- "Multi-class" 
      pred.ord <- paste(levels(as.factor(predictor)), collapse = ", ")
    }
  }
  if(is.null(paired)) pair.det <- "No" else pair.det <- "Yes"
  
  det <- data.frame(Features = nrow(count_table),
                               Samples = ncol(count_table),
                               Predictor = pred.det,
                               Ordering = pred.ord,
                               Paired = pair.det,
                               Covars = paste(names(covars), collapse = ", "),
                               Relative = relative,
                               OutAll = out.all)
  rownames(det) <- ""
  output.details <- as.data.frame(t(det))
  colnames(det) <- ""
  
  return(list(raw = df.raw, adj = df.adj, est = df.est,details = det, results = newresults))

}




