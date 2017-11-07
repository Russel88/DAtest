#' Plot association between abundance of a feature and predictor
#'
#' Plot association between abundance of a feature and predictor, modified if paired and covars are available
#' 
#' Boxplot for categorical variables, points and smooth line for quantitative variable.
#' If a paired variable is supplied, it is always plotted as points with lines grouped by the paired variable
#' If covars are supplied data is split in facets. Quantitative covars are cut in intervals according to the quantiles given in covar.quant
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2", "zig" and "fri"
#' @param covars Either a named list with covariates, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param feature Name of feature to plot. Should be in rownames of data (or taxa_names(data) if data is a phyloseq object)
#' @param relative Logical. If TRUE (default) abundances are made relative
#' @param logScale Logical. Should abundances be log10-scaled? After normalization if relative is TRUE
#' @param delta Pseudocount for log10 normalization
#' @param covar.quant Quantiles for cutting quantitative covars
#' @return A ggplot
#' @export

featurePlot <- function(data, predictor, paired = NULL, covars = NULL, feature = NULL, relative = TRUE, logScale = FALSE, delta = 0.001, covar.quant = c(0,1/3,2/3,1)){

  stopifnot(exists("data"),exists("predictor"))

  library(ggplot2)
  
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

  if(is.null(feature)) stop("feature cannot be NULL, it should match a rownames in data")
  if(!feature %in% rownames(count_table)) stop("feature not found in data")
  
  # Relative and logscale
  if(relative) count_table <- apply(count_table, 2, function(x) x/sum(x))
  if(logScale) count_table <- log10(count_table + delta)
  
  # Dataframe
  df <- data.frame(Abundance = count_table[rownames(count_table) == feature,],
                   Predictor = predictor)
  if(!is.null(paired)){
    df <- cbind(df, paired)
    colnames(df) <- c("Abundance","Predictor","Paired")
  }
  if(!is.null(covars)){
    oldcn <- colnames(df)
    for(i in seq_along(covars)){
      subco <- covars[[i]]
      if(is.numeric(subco) & length(unique(subco)) > 3){
        subco.new <- cut(subco, breaks = as.numeric(quantile(subco, probs = covar.quant)), include.lowest = TRUE)
        df <- cbind(df, subco.new)
      } else {
        df <- cbind(df, subco)
      }
    }
    colnames(df) <- c(oldcn,names(covars))
  }
  
  # Covar facetting
  if(!is.null(covars)){
    cv.fac <- as.formula(paste("~",paste(names(covars),collapse = "+")))
  }
  
  # Plotting
  if(is.null(covars)){
    if(is.null(paired)){
      p <- ggplot(df, aes(Predictor, Abundance)) +
        theme_bw()
      if(is.numeric(predictor)) {p <- p + geom_point() + geom_smooth()} else {p <- p + geom_boxplot()}
    } else {
      if(all(table(df$Predictor,df$Paired)==1)){
        p <- ggplot(df, aes(Predictor, Abundance, group = Paired, colour = Paired)) +
          theme_bw() +
          geom_line() 
      } else {
        p <- ggplot(df, aes(Predictor, Abundance, group = Paired, colour = Paired)) +
          theme_bw() +
          geom_point() +
          stat_summary(fun.y=mean, geom="line") 
      }
    }
  } else {
    if(is.null(paired)){
      p <- ggplot(df, aes(Predictor, Abundance)) +
        theme_bw() +
        facet_wrap(cv.fac)
      if(is.numeric(predictor)) {p <- p + geom_point() + geom_smooth()} else {p <- p + geom_boxplot()}
    } else {
      if(all(table(df$Predictor,df$Paired)==1)){
        p <- ggplot(df, aes(Predictor, Abundance, group = Paired, colour = Paired)) +
          theme_bw() +
          geom_line() +
          facet_wrap(cv.fac)
      } else {
        p <- ggplot(df, aes(Predictor, Abundance, group = Paired, colour = Paired)) +
          theme_bw() +
          geom_point() +
          stat_summary(fun.y=mean, geom="line") +
          facet_wrap(cv.fac)
      }
    }
  }
  
  # Title
  if(class(data) == "phyloseq"){
    tax <- tax_table(data)
    subtax <- tax[rownames(tax) == feature,]
    p <- p + ggtitle(feature, subtitle = paste(subtax, collapse = "_"))
  } else {
    p <- p + ggtitle(feature)
  }

  p
}




