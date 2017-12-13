#' Are some groups/taxa overrepresented among significant features
#'
#' Test if some groups of features are overpresented among significant features. 
#' The groups can be anything; for OTU data e.g. genera/family/order/class/phylum, for transciptomics e.g. KEGG pathway.
#' 
#' OR in output is odds ratio from fisher's exact test. If OR is above 1 it means that the group is overrepresented among significant features.
#' @param results Data.frame with results from a \code{DAtest} function
#' @param group.df Data.frame with columns defining the groups. rownames should name the features matching the \code{Feature} column in \code{results}. E.g. \code{tax_table} from a \code{phyloseq} object
#' @param group.cols Numeric vector defining which column(s) contain(s) the groups in \code{group.df}. Default first column.
#' @param split If TRUE will split tests in positive and negative effect sizes if possible. Default TRUE
#' @param alpha Threshold for significance calling. Default 0.05
#' @param p.adj Method for p-value adjustment. Default "fdr"
#' @param alternative What to test for. "greater" (default) is testing only overrepresentation (OR > 1), "less" only underrepresentation (OR < 1), and "two.sided" tests over- and under-representation (OR != 1)
#' @return A data.frame with odds ratios (OR), p-values, adjusted p-values, groups, name of groups, and direction of effect if split = TRUE
#' @export
groupSig <- function(results, group.df, group.cols = 1, split = TRUE, alpha = 0.05, p.adj = "fdr", alternative = "greater"){

  stopifnot(exists("results"),exists("group.df"))

  # Name of estimate for different methods
  est.name <- list(sam = c("ordering","log2FC"),
                   znb = c("ordering","log2FC"),
                   zpo = c("ordering","log2FC"),
                   qpo = c("ordering","log2FC"),
                   poi = c("ordering","log2FC"),
                   neb = c("ordering","log2FC"),
                   lrm = c("ordering","log2FC"),
                   llm = c("ordering","log2FC"),
                   llm2 = c("ordering","log2FC"),
                   vli = c("ordering","logFC"),
                   lim = c("ordering","logFC"),
                   lli = c("ordering","logFC"),
                   lli2 = c("ordering","logFC"),
                   pea = "cor",
                   spe = "rho",
                   per = c("ordering","log2FC"),
                   bay = "ordering",
                   adx.t = c("ordering","effect"),
                   adx.w = c("ordering","effect"),
                   wil = c("ordering","log2FC"),
                   ttt = c("ordering","log2FC"),
                   ltt = c("ordering","log2FC"),
                   ltt2 = c("ordering","log2FC"),
                   erq = c("ordering","log2FC"),
                   ere = c("ordering","log2FC"),
                   erq2 = c("ordering","log2FC"),
                   ere2 = c("ordering","log2FC"),
                   msf = c("ordering","log2FC"),
                   zig = c("ordering","log2FC"),
                   ds2 = c("ordering","log2FoldChange"),
                   ds2x = c("ordering","log2FoldChange"),
                   rai = c("ordering","log2FC"))
  
  # The method
  method <- unique(gsub("\\)","",gsub(".*\\(","",results$Method)))
  
  # Subset group dataframe
  group.df <- as.data.frame(group.df)[,group.cols, drop = FALSE]
  
  # Fix for anc
  if(method == "anc"){
    results$pval.adj <- 1
    results[results$Detected == "Yes","pval.adj"] <- 0
  }
  
  # Fix for sam
  if(method == "sam"){
    if("Sig" %in% colnames(results)){
      results$pval.adj <- 1
      results[results$Sig == "Yes","pval.adj"] <- 0
    } else {
      results$pval.adj <- 1
      results[results$Sig.up == "Yes","pval.adj"] <- 0
      results[results$Sig.lo == "Yes","pval.adj"] <- 0
    }
  }
  
  # Significant
  feat.sig <- na.omit(results[results$pval.adj < alpha,"Feature"])
  if(length(feat.sig) == 0) stop("No significant features")
  group.df$Sig <- 0
  group.df[rownames(group.df) %in% feat.sig,"Sig"] <- 1
  if(all(group.df$Sig == 1)) stop("All features are significant!")
  
  # The effect size name
  split.name <- tryCatch(est.name[[method]], error = function(e) NULL)
  if(is.null(split.name)) split.name <- "DAerror"
  if(length(split.name) == 2){
    if("ordering" %in% colnames(results)) split.name <- split.name[1] else split.name <- split.name[2]
  } 
  
  # For splitting in positive and negative
  if(method %in% names(est.name) & split.name %in% colnames(results) & split){
    
    # Find positive and negatives
    if(split.name != "ordering"){
      pos.feat <- results[results[split.name] > 0,"Feature"]
      neg.feat <- results[results[split.name] < 0,"Feature"]
    } else {
      pos.feat <- results[results$ordering == unique(results$ordering)[1],"Feature"]
      neg.feat <- results[results$ordering == unique(results$ordering)[2],"Feature"]
      if(length(unique(results$ordering)) == 3) other.feat <- na.omit(results[results$ordering == unique(results$ordering)[3],"Feature"])
    }
    pos.feat <- na.omit(pos.feat)
    neg.feat <- na.omit(neg.feat)
   
    # Put effect size in data.frame
    group.df$Effect <- NA
    if(split.name != "ordering"){
      group.df[rownames(group.df) %in% pos.feat,"Effect"] <- "Positive"
      group.df[rownames(group.df) %in% neg.feat,"Effect"] <- "Negative"
    } else {
      group.df[rownames(group.df) %in% pos.feat,"Effect"] <- as.character(unique(results$ordering)[1])
      group.df[rownames(group.df) %in% neg.feat,"Effect"] <- as.character(unique(results$ordering)[2])
      if(length(unique(results$ordering)) == 3) group.df[rownames(group.df) %in% other.feat,"Effect"] <- as.character(unique(results$ordering)[3])
    }
    group.df <- na.omit(group.df)

    # Run the tests
    allres <- list()
    for(j in 1:(ncol(group.df)-2)){
      subres <- data.frame(Group = NA,
                           OR = NA,
                           pval = NA)
      group.df[,j] <- paste0(group.df[,j],"_DAtest_",group.df$Effect)
      # Table
      tab <- table(group.df[,j],group.df$Sig)
      for(i in 1:nrow(tab)){
        # The test
        fit <- fisher.test(rbind(colSums(tab[-i,]),tab[i,]), alternative = alternative)
        subres[i,"Group"] <- rownames(tab[i,,drop = FALSE])
        subres[i,"OR"] <- fit$estimate
        subres[i,"pval"] <- fit$p.value
      }
      subres$GroupName <- colnames(group.df)[j]
      subres$Effect <- gsub(".*_DAtest_","",subres$Group)
      subres$Group <- gsub("_DAtest_.*","",subres$Group)
      allres[[j]] <- subres
    }
    final <- do.call(rbind,allres)
    final <- final[,c(4,1,5,2,3)]
    
  } else {
    message("Note: Results not splitted in positive and negativ effects")
    allres <- list()
    for(j in 1:(ncol(group.df)-1)){
      subres <- data.frame(Group = NA,
                           OR = NA,
                           pval = NA)
      # Table
      tab <- table(group.df[,j],group.df$Sig)
      for(i in 1:nrow(tab)){
        # Collapse table
        fit <- fisher.test(rbind(colSums(tab[-i,]),tab[i,]), alternative = alternative)
        subres[i,"Group"] <- rownames(tab[i,,drop = FALSE])
        subres[i,"OR"] <- fit$estimate
        subres[i,"pval"] <- fit$p.value
      }
      subres$GroupName <- colnames(group.df)[j]
      allres[[j]] <- subres
    }
    final <- do.call(rbind,allres)
    final <- final[,c(4,1,2,3)]
  }
  
  final$pval.adj <- p.adjust(final$pval, method = p.adj)
  return(final)
}






