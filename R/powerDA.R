#' Estimating (empirical) statistical power
#'
#' Estimating (empirical) statistical power for a specific differential abundance and expression method on a specific dataset
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. Only for "anc", "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "ds2x", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param test Character. Which test to include. See \code{testDA} for details on the implemented tests. 
#' @param effectSizes Numierc. The effect sizes for the spike-ins. Default \code{c(2,4,8,16,32)}
#' @param alpha.p p-value threshold for false positive rates. Default 0.05
#' @param alpha.q q-value threshold for determining significance for \code{empirical power}. Default 0.05. This will change \code{fdr.output} for "sam" and \code{sig} for "anc". 
#' @param p.adj Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"
#' @param R Integer. Number of times to run the tests. Default 5
#' @param relative Logical. If TRUE (default) abundances are made relative for "ttt", "ltt", "wil", "per", "aov", "lao", "kru", "lim", "lli", "lrm", "llm", "spe" and "pea", and there is an offset of \code{log(LibrarySize)} for "neb", "poi", "qpo", "zpo" and "znb"
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. \code{k=c(5,10,15)}: 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default NULL, which will spike 1 percent of the total amount of features in each tertile (a total of 3 percent)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param rng.seed Numeric. Seed for reproducibility. Default 123
#' @param args List. A list with lists of arguments passed to the different methods. See details for more.
#' @param out.all If TRUE linear models will output results and p-values from \code{anova}/\code{drop1}, ds2/ds2x will run LRT and not Wald test, erq and erq2 will produce one p-value for the predictor, and lim, lli, lli2, lim, vli will run F-tests. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param core.check If TRUE will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.  
#' @param verbose If TRUE will print informative messages
#' @details Currently implemented methods: see \code{testDA}
#' @return An object of class \code{DAPower}, which contains a list with 1: A data.frame with results, 2: alpha.p value, 3: alpha.q values
#' @import snow doSNOW foreach utils
#' @importFrom parallel detectCores
#' @importFrom pROC roc
#' @export

powerDA <- function(data, predictor, paired = NULL, covars = NULL, test = NULL, effectSizes = c(2,4,8,16,32), alpha.p = 0.05, alpha.q = 0.05, p.adj = "fdr", R = 5, relative = TRUE, k = NULL, cores = (detectCores()-1), rng.seed = 123, args = list(), out.all = NULL, core.check = TRUE, verbose = TRUE){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(core.check){
    if(cores > 20){
      ANSWER <- readline(paste("You are about to run testDA using",cores,"cores. Enter y to proceed "))
      if(ANSWER != "y") stop("Process aborted")
    }
  }

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
  if(is.null(test)) stop("'test' has to be specified")
  if(relative) if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) stop("count_table must only contain integer values")
  if(min(count_table) < 0) stop("count_table contains negative values")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(unique(predictor)) < 2) stop("predictor should have at least two levels")
  if(length(test) != 1) stop("'test' has to have length 1")
  
  # Set seed
  set.seed(rng.seed)
  if(verbose) message(paste("Seed is set to",rng.seed))

  # Create some random seeds for each run
  seeds <- rpois(R, lambda = (1:R)*1e6)

  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Spike vs no features
  if(is.null(k)){
    k <- rep(round(nrow(count_table)*0.01),3)
  } 
  if(sum(k) == nrow(count_table)) stop("Set to spike all features. Change k argument")
  if(sum(k) > nrow(count_table)) stop("Set to spike more features than are present in the data. Change k argument")
  if(sum(k) == 0) k <- c(1,1,1)                                  
  if(sum(k) < 15 & sum(k) >= 10 & R <= 10) message("Few features spiked. Increase 'k' or set 'R' to more than 10 to ensure proper estimations")
  if(sum(k) < 10 & sum(k) >= 5 & R <= 20) message("Few features spiked. Increase 'k' or set 'R' to more than 20 to ensure proper estimations")                                  
  if(sum(k) < 5 & R <= 50) message("Very few features spiked. Increase 'k' or set 'R' to more than 50 to ensure proper estimations")
                                    
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
  spikeds.l <- list()
  for(eff in seq_along(effectSizes)){
    spikeds.l[[eff]] <- lapply(1:R,function(x) spikein(count_table, rands[[x]], effectSizes[eff],  k, num.pred, relative))
  }
  spikeds <- do.call(c, spikeds.l)
  count_tables <- lapply(1:(R*length(effectSizes)),function(x) spikeds[[x]][[1]])
  
  # Test list
  tests.par <- paste0(unlist(lapply(effectSizes,function(x) rep(x,R))),"-",rep(paste0(1:R,"_",rep(test,R)),R))

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
  results <- foreach(i = tests.par, .options.snow = opts) %dopar% {

    # Extract run info
    what.run <- which(i == tests.par)
    run.no <- as.numeric(gsub(".*-","",gsub("_.*","",i)))
    i <- gsub(".*_","",i)

    # Set subseed
    set.seed(seeds[run.no])
    
    # Extract test arguments
    if(!all(names(args) %in% test)) stop("One or more names in list with additional arguments does not match name of test")
    for(j in seq_along(args)){
      assign(paste0(names(args)[j],".DAargs"),args[[j]],pos=1)
    }
    test.args <- paste0(test,".DAargs")
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
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars),zzz.DAargs)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj),wil.DAargs)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj),ttt.DAargs)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj),ltt.DAargs)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,p.adj),ltt2.DAargs)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj),neb.DAargs)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),erq.DAargs)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj),ere.DAargs)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),erq2.DAargs)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj),ere2.DAargs)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj),msf.DAargs)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,p.adj),zig.DAargs)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),ds2.DAargs)),
                               ds2x = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),ds2x.DAargs)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj),per.DAargs)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,p.adj),bay.DAargs)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]]),adx.DAargs)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj),lim.DAargs)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj),lli.DAargs)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),lli2.DAargs)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj),kru.DAargs)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,p.adj),aov.DAargs)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,p.adj),lao.DAargs)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,p.adj),lao2.DAargs)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars, relative,out.all,p.adj),lrm.DAargs)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj),llm.DAargs)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),llm2.DAargs)),
                               rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj),rai.DAargs)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj),spe.DAargs)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj),pea.DAargs)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj),poi.DAargs)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj),qpo.DAargs)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj),vli.DAargs)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj),zpo.DAargs)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj),znb.DAargs)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj),fri.DAargs)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj),qua.DAargs)),
                               anc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,sig = alpha.q),anc.DAargs)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,fdr.output = alpha.q),sam.DAargs))),
                        
                        error = function(e) NULL)
    
    if(!i %in% c("anc","sam","adx")){
      res.sub[is.na(res.sub$pval),"pval"] <- 1
      res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
    }
    
    return(res.sub)
    
  }
  names(results) <- tests.par
  
  # Split aldex
  if(test == "adx"){
    results.t <- lapply(results, function(x) x[,c(1:9,12,13)])
    results.w <- lapply(results, function(x) x[,c(1:7,10:13)])
    for(k in seq_along(results.t)){
      colnames(results.t[[k]]) <- c(colnames(results[[1]])[1:7],"pval","pval.adj","ordering","Feature")
      colnames(results.w[[k]]) <- c(colnames(results[[1]])[1:7],"pval","pval.adj","ordering","Feature")
      results.t[[k]]$Method <- "ALDEx2 t-test (adx)"
      results.w[[k]]$Method <- "ALDEx2 wilcox (adx)"
    }
    names(results.t) <- paste0(names(results.t),".t")
    names(results.w) <- paste0(names(results.w),".w")
    results <- c(results.t,results.w)
    tests.par <- names(results)
    spikeds <- c(spikeds,spikeds)
  }
  
  final.results <- foreach(r = tests.par, .combine = rbind) %do% {

    res.sub <- results[names(results) == r][[1]]
    
    # Make pseudo-pval for ANCOM
    if(test == "anc"){
      suppressWarnings(min.d <- min(res.sub[res.sub$Detected == "Yes","W"]))
      if(min.d == Inf) min.d <- max(res.sub$W)+1
      res.sub$pval <- 1/(res.sub$W+1) * alpha.p/(1/(min.d+1))
      res.sub$pval.adj <- 1/(res.sub$W+1) * alpha.q/(1/(min.d+1))
    }
    
    # Make pseudo-pval for SAMseq
    if(test == "sam"){
      res.sub$pval <- 1/rank(res.sub$Score)
      res.sub$pval.adj <- 1
      res.sub[res.sub$Sig == "Yes","pval.adj"] <- 0 
    }
    
    # Confusion matrix
    totalPos <- nrow(res.sub[res.sub$pval <= alpha.p,])
    totalNeg <- nrow(res.sub[res.sub$pval > alpha.p,])
    truePos <- sum(res.sub[res.sub$pval <= alpha.p,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falseNeg <- sum(res.sub[res.sub$pval > alpha.p,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR
    if(test == "sam"){
      fpr <- NA
    } else {
      if((falsePos + trueNeg) != 0){
        fpr <- falsePos / (falsePos + trueNeg)
      } else {
        fpr <- 0
      }
    }

    # Confusion matrix adjusted
    totalPos.adj <- nrow(res.sub[res.sub$pval.adj <= alpha.q,])
    truePos.adj <- sum(res.sub[res.sub$pval.adj <= alpha.q,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falsePos.adj <- totalPos.adj - truePos.adj
    
    # FDR 
    if(totalPos.adj != 0){
      fdr <- falsePos.adj / totalPos.adj
    } else {
      fdr <- 0
    }
    
    # Spike detection rate (empircal power aka sensitivity)
    sdr <- truePos.adj / sum(k)
    
    # AUC
    test_roc <- NULL
     tryCatch(
       test_roc <- pROC::roc(as.numeric(res.sub$Feature %in% spikeds[[which(r == tests.par)]][[2]]) ~ res.sub$pval, auc=TRUE, direction = ">"),
       error = function(e) NULL)
     if(!is.null(test_roc)){
       auc <- as.numeric(test_roc$auc) 
     } else {
       auc <- 0.5
     }

    # Combine and return
    df.combined <- data.frame(Method = res.sub$Method[1],
                              Run = as.numeric(gsub(".*-","",gsub("_.*","",r))),
                              EffectSize = as.numeric(gsub("-.*","",gsub("_.*","",r))),
                              AUC = auc,
                              FPR = fpr,
                              FDR = fdr,
                              SDR = sdr)
    rownames(df.combined) <- NULL

    return(df.combined)
    
  }
  final <- list(final.results,alpha.p,alpha.q)
  class(final) <- "DAPower"
  return(final)
}
