#' Comparing differential abundance/expression methods by Empirical power and False Discovery Rate
#'
#' Calculating Power, False Discovery Rates, False Positive Rates and AUC (Area Under the Receiver Operating Characteristic (ROC) Curve) for various differential abundance and expression methods
#' 
#' mva is excluded by default, as it is slow.
#' @param data Either a data.frame with counts/abundances, OR a \code{phyloseq} object. If a data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation.
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param R Integer. Number of times to run the tests. Default 20
#' @param tests Character. Which tests to include. Default all
#' @param relative Logical. TRUE (default) for compositional data. FALSE for absolute abundances or pre-normalized data.
#' @param effectSize Numeric. The effect size for the spike-ins. Default 5
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. \code{k=c(5,10,15)}: 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default NULL, which will spike 2 percent of the total amount of features in each tertile (a total of 6 percent), but minimum c(5,5,5)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param p.adj Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.all If TRUE linear models will output results and p-values from \code{anova}/\code{drop1}, ds2/ds2x will run LRT and not Wald test, erq and erq2 will produce one p-value for the predictor, and limma will run F-tests. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param alpha q-value threshold for determining significance for \code{Power}. Default 0.1
#' @param core.check If TRUE will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.  
#' @param verbose If TRUE will print informative messages
#' @return An object of class \code{DA}, which contains a list of results:
#' \itemize{
#'  \item table - FPR, AUC and spike detection rate for each run
#'  \item results - A complete list of output from all the methods. Example: Get wilcoxon results from 2. run as such: \code{$results[[2]]["wil"]}
#'  \item details - A dataframe with details from the run
#'  \item run.times - A dataframe with average run times of the different methods
#' }
#' 
#' @import stats snow doSNOW foreach utils doParallel
#' @importFrom parallel detectCores
#' @importFrom pROC roc
#' @export

testDA <- function(data, predictor, paired = NULL, covars = NULL, R = 20,
                   tests = c("bay","ds2","ds2x","per","adx","znb","zpo","msf","zig",
                             "erq","erq2","neb","qpo","poi","sam",
                             "lrm","llm","llm2","lma","lmc",
                             "ere","ere2","pea","spe",
                             "wil","kru","qua","fri",
                             "ttt","ltt","ltt2","tta","ttc",
                             "aov","lao","lao2","aoa","aoc",
                             "vli","lim","lli","lli2","lia","lic"),
                   relative = TRUE, effectSize = 5, k = NULL, cores = (detectCores()-1),
                   p.adj = "fdr", args = list(), out.all = NULL, alpha = 0.1, core.check = TRUE, verbose = TRUE){

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
    for(i in seq_along(covars)){
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
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0 && verbose) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  if(nrow(count_table) <= 15) warning("Dataset contains very few features") 
  
  # Prune tests argument
  decimal <- zeroes <- FALSE
  if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) decimal <- TRUE
  if(any(count_table == 0)) zeroes <- TRUE
  tests <- unique(tests)
  if(!"zzz" %in% tests) tests <- prune.tests.DA(tests, predictor, paired, covars, relative, decimal, zeroes)
  tests.par <- paste0(unlist(lapply(seq_len(R), function(x) rep(x,length(tests)))),"_",rep(tests,R))
  if(length(tests) == 0) stop("No tests to run!")
  
  # Run time warnings
  if(verbose){
    if("neb" %in% tests & !is.null(paired)){
      message("As 'neb' is included and a 'paired' variable is supplied, this might take a long time")
    }
  }

  if(verbose) message(paste("Running on",cores,"cores"))

  # Spike vs no features
  if(is.null(k)){
    k <- rep(round(nrow(count_table)*0.02),3)
    if(sum(k) < 15){
      k <- c(5,5,5)
    } 
  } 
  if(sum(k) == nrow(count_table)) stop("Set to spike all features. Can't calculate FPR or AUC. Change k argument")
  if(sum(k) > nrow(count_table)) stop("Set to spike more features than are present in the data. Change k argument")
  if(sum(k) < 15 & sum(k) >= 10 & R <= 10) warning("Few features spiked. Increase 'k' or set 'R' to more than 10 to ensure proper estimation of AUC and FPR")
  if(sum(k) < 10 & sum(k) >= 5 & R <= 20) warning("Few features spiked. Increase 'k' or set 'R' to more than 20 to ensure proper estimation of AUC and FPR")
  if(sum(k) < 5 & R <= 50) warning("Very few features spiked. Increase 'k' or set 'R' to more than 50 to ensure proper estimation of AUC and FPR")
  if(sum(k) > nrow(count_table)/2) warning("Set to spike more than half of the dataset, which might give unreliable estimates, Change k argument")
                                 
  # predictor
  if(any(is.na(predictor))) warning("Predictor contains NAs!")
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
  if(!is.null(paired)){
    if(verbose) message(paste("The paired variable has",length(unique(paired)),"levels"))
    if(length(unique(paired)) < 5) warning("The paired variable has less than 5 levels. Mixed-effect models are excluded")
  }
  
  # out.all
  if(is.null(out.all)){
    if(length(unique(predictor)) == 2) out.all <- FALSE
    if(length(unique(predictor)) > 2) out.all <- TRUE
    if(num.pred) out.all <- FALSE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in seq_along(covars)){
      if(verbose) if(any(is.na(covars[[i]]))) warning(names(covars)[i],"contains NAs!")
      if(is.numeric(covars[[i]][1])){
        if(verbose) message(paste(names(covars)[i],"is assumed to be a quantitative variable, ranging from",min(covars[[i]], na.rm = TRUE),"to",max(covars[[i]], na.rm = TRUE)))
      } else {
        if(verbose) message(paste(names(covars)[i],"is assumed to be a categorical variable with",length(unique(covars[[i]])),"levels:",paste(levels(as.factor(covars[[i]])),collapse = ", ")))
      }
    }
  }
  
  if(verbose) cat("Spikeing...\n")
  # Shuffle predictor
  if(is.null(paired)){
    rands <- lapply(seq_len(R),function(x) sample(predictor))
  } else {
    rands <- lapply(seq_len(R),function(x) unsplit(lapply(split(predictor,paired), sample), paired))
  }
  
  # Spikeins
  spikeds <- lapply(seq_len(R),function(x) spikein(count_table, rands[[x]], effectSize,  k, num.pred, relative))
  count_tables <- lapply(seq_len(R),function(x) spikeds[[x]][[1]])
  
  # Extract test arguments
  if(!all(names(args) %in% tests)) stop("One or more names in list with additional arguments does not match names of tests")
  argsL <- list()
  for(a in seq_along(tests)){
    argsL[tests[a]] <- args[tests[a]]
  }
  
  ### Run tests
  if(verbose) cat(paste("Testing", length(tests),"methods", R,"times each...\n"))
  # Progress bar
  pb <- txtProgressBar(max = length(tests.par), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Start parallel
  if(cores == 1) {
    registerDoSEQ() 
  } else {
    cl <- parallel::makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  }

  # Run the tests in parallel
  results <- foreach(i = tests.par, .options.snow = opts) %dopar% {
    
    t1.sub <- proc.time()
    
    # Extract run info
    run.no <- as.numeric(gsub("_.*","",i))
    i <- gsub(".*_","",i)

    # zzz tests
    if(!is.na(pmatch("zzz",i))){
      j <- i
      i <- "zzz"
    }
    
    # Run tests
    res.sub <- tryCatch(switch(i,
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars), argsL[[j]])),
                               mva = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative, p.adj), argsL[[i]])),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj), argsL[[i]])),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj), argsL[[i]])),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj), argsL[[i]])),
                               tta = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, p.adj), argsL[[i]])),
                               ttc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, p.adj), argsL[[i]])),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, p.adj), argsL[[i]])),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj), argsL[[i]])),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj), argsL[[i]])),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj), argsL[[i]])),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], p.adj), argsL[[i]])),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars, p.adj), argsL[[i]])),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               ds2x = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired, relative, p.adj), argsL[[i]])),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]), argsL[[i]])),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]]), argsL[[i]])),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj), argsL[[i]])),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj), argsL[[i]])),
                               lia = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               lic = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]], relative, p.adj), argsL[[i]])),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, relative, p.adj), argsL[[i]])),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative, p.adj), argsL[[i]])),
                               aoa = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, p.adj), argsL[[i]])),
                               aoc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, p.adj), argsL[[i]])),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars, p.adj), argsL[[i]])),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars, relative,out.all, p.adj), argsL[[i]])),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj), argsL[[i]])),
                               lma = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               lmc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative, p.adj), argsL[[i]])),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],relative, p.adj), argsL[[i]])),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,relative,out.all, p.adj), argsL[[i]])),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj), argsL[[i]])),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,covars,out.all, p.adj), argsL[[i]])),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj), argsL[[i]])),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],covars,relative,out.all, p.adj), argsL[[i]])),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj), argsL[[i]])),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,relative, p.adj), argsL[[i]])),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[run.no]],rands[[run.no]],paired,fdr.output = alpha), argsL[[i]]))),
                        error = function(e) NULL)
    
    if(!is.null(res.sub) & (!i %in% c("sam","adx"))){
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
  results <- results[!vapply(results,is.null)]
  
  cat("\n")
  if(length(unique(gsub(".*_","",names(results)))) != length(tests)){
    if(length(tests) - length(unique(gsub(".*_","",names(results)))) == 1){
      if(verbose) message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"was excluded due to failure"))
    } else {
      if(verbose) message(paste(paste(tests[!tests %in% unique(gsub(".*_","",names(results)))],collapse = ", "),"were excluded due to failure"))
    }

    # Produce informative messages
    if(all(tests[!tests %in% unique(gsub(".*_","",names(results)))] == "sam")){
      if(verbose) message("sam usually fails if some samples has too many zeroes")
    }
    if(all(c("sam","ere2","erq2","ds2x") %in% tests[!tests %in% unique(gsub(".*_","",names(results)))])){
      if(verbose) message("sam, ere2, erq2 and ds2x usually fails if all features contain at least one zero")
    }
      
    tests <- unique(gsub(".*_","",names(results)))
  }
  
  r <- NULL
  final.results <- foreach(r = seq_len(R)) %do% {

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
    rsp <- NULL
    res.sub <- foreach(rsp = seq_along(res.sub)) %do% {
      temp <- res.sub[[rsp]]
      temp$Spiked <- "No"
      temp[temp$Feature %in% spikeds[[r]][[2]],"Spiked"] <- "Yes"
      return(temp)
    }
    names(res.sub) <- newnames
    
    # Confusion matrix
    totalPos <- vapply(res.sub,function(x) nrow(x[x$pval <= 0.05,]))
    totalNeg <- vapply(res.sub,function(x) nrow(x[x$pval > 0.05,])) 
    trueNeg <- totalNeg  #if effectSize == 1
    truePos <- 0  #if effectSize == 1
    falseNeg <- 0 #if effectSize == 1
    if(effectSize != 1){
      truePos <- vapply(res.sub, function(x) sum(x[x$pval <= 0.05,"Feature"] %in% spikeds[[r]][[2]]))
      falseNeg <- vapply(res.sub, function(x) sum(x[x$pval > 0.05,"Feature"] %in% spikeds[[r]][[2]]))
    }
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR 
    fprs <- vapply(seq_along(res.sub), function(x) {
      if((falsePos[x] + trueNeg[x]) != 0){
        falsePos[x] / (falsePos[x] + trueNeg[x])
      } else {0}})
    
    
    # Spike detection rate
    # True positive for adjusted p-values
    truePos.adj <- 0  #if effectSize == 1
    if(effectSize != 1){
      truePos.adj <- vapply(res.sub, function(x) sum(x[x$pval.adj <= alpha,"Feature"] %in% spikeds[[r]][[2]]))
    }
    sdrs <- vapply(seq_along(res.sub), function(x) truePos.adj[x] / sum(k))
    
    # AUC
    aucs <- vapply(seq_along(res.sub), function(x) {
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
    
    # Confusion matrix adjusted
    totalPos.adj <- vapply(res.sub, function(x) nrow(x[x$pval.adj <= alpha,]))
    truePos.adj <- vapply(res.sub, function(x) sum(x[x$pval.adj <= alpha,"Feature"] %in% spikeds[[r]][[2]]))
    falsePos.adj <- totalPos.adj - truePos.adj
    
    fdrs <- vapply(seq_along(res.sub), function(x) {
      if(totalPos.adj[x] != 0){
        falsePos.adj[x] / totalPos.adj[x]
      } else {0}})
    
    # Combine and return
    df.combined <- data.frame(Method = vapply(res.sub, function(x) x$Method[1]),
                              AUC = aucs,
                              FPR = fprs,
                              FDR = fdrs,
                              Power = sdrs,
                              Run = r)
    rownames(df.combined) <- NULL

    # SAMseq FPRs not possible
    if("sam" %in% newnames){
      df.combined[df.combined$Method == "SAMseq (sam)","FPR"] <- NA
    }
        return(list(df.combined,res.sub))
  }
  
  output.results <- do.call(rbind,lapply(final.results, function(x) x[[1]]))
  output.all.results <- lapply(final.results, function(x) x[[2]])
  
  if(num.pred){
    pred.det <- "Quantitative"
    pred.ord <- paste(min(predictor,na.rm=TRUE),"to",max(predictor,na.rm=TRUE))
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
  
  return(out)
}
