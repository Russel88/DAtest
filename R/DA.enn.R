#' ENNB from https://cals.arizona.edu/~anling/
#' @import glmnet MASS
#' @export

DA.enn <- function(count_table, outcome, TMM.option = 1, p.adj){
  
  #################################################################################
  #                                                                               #
  #  ENNB: A two-stage statistical procedure for feature selection and comparison 
  #  in fuctional analysis of metagenomes                                            #
  #				                                         	#
  #  input:									#
  #   X: Matrix of count data 							#
  #   Y: Phenotype info of samples                                                #
  #  fileoutput: the outout file name						#
  #   TMM.option: 1 or 2. option of "1" is for an approach using the mean of the	#
  #   effective library sizes as a reference library size in TMM normalization; 	#
  #   option of "2" represents an approach to regenerating counts 		#
  #       with a common dispersion. 						#
  #                                                                               #
  # output:                                                                       #                    
  # summary table of the final results after two stages.                   	#                                                             
  #################################################################################
  
  TwoStage_Package <- function(X,Y, TMM.option){
    require(glmnet)   # Elastic-net
    require(MASS)     # GLMs-NB
    require(edgeR)    # TMM normalization
    
    threshold = 0.05  # FDR to control the type I error at significance level of 0.05
    
    #####################################################
    ###   Trimmed Mean (TMM) (Based on edgeR Paper)   ###
    #####################################################
    TMMNorm = function(X, Y, TMM.option){
      factors = calcNormFactors(X,method="TMM") #Calculate normaization factors
      if(TMM.option==1){
        eff.lib.size = colSums(X)*factors;
        ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size 
        X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts
      } else if(TMM.option==2){ #Regenerate counts with a common dispersion.
        group = as.factor(Y);
        X.tmp = DGEList(counts=X, group=group);
        X.tmp = calcNormFactors(X.tmp);
        X.tmp = estimateCommonDisp(X.tmp);
        X.output = round(as.matrix(X.tmp$pseudo.counts));
      } else {
        stop("TMM.option must be either 1 or 2!");
      }
      return(X.output);
    }
    
    ############################################################################
    #### Negative binomial generalized linear model for two group comparison ###
    ############################################################################
    get.glmNB <- function(count,category){
      p.val = rep(NA,ncol(count))  #Declare initial pvalue
      
      for (i in 1:ncol(count)) {
        dat = data.frame(testing=count[,i], groups=as.factor(category))
        if (sum(dat$testing) != 0) { #in case that all counts are zero
          fit = glm.nb(testing ~ groups, data=dat)
          p.val[i] = summary(fit)$coefficients[2,4] 
        }
      }
      p.val[is.na(p.val)] = 1
      p.adj = p.adjust(p.val, method="fdr") #Adjusted p-value 
      res = cbind(p.val,p.adj)
      rownames(res) = colnames(count)
      return(res)  
    }
    
    ##############################################################################
    ###                           Two-Stage Procedure                          ###
    ##############################################################################
    #################################
    ###   Two group comparison   ####
    #################################
    Twostage <- function(X, Y){
      ### A matrix contains cvm error information
      cv.mat = matrix(NA, nrow = length(alpha), ncol = 3)
      colnames(cv.mat) = c("Alpha", "Lambda", "CVM")
      cv.mat[,1] = alpha
      
      #### Cross-validation to determine the optimal alpha and lambda giving minimum cross-validation error 
      npop1 = table(Y)[1]
      npop2 = table(Y)[2]
      N = npop1 + npop2
      if (N <= 15) nfold = 3
      if (N > 15  & N < 50) nfold = 5
      if (N >= 50) nfold = 10  
      
      ### Fit a GLM with elastic net regularization
      for (iter in 1:3){
        set.seed(iter)
        
        #Select a random number of subjects
        newdata.ind = sample(1:nrow(X), floor(0.9*nrow(X)))
        
        #Subset the data
        X.new = X[newdata.ind,]
        Y.new = Y[newdata.ind]
        
        # Cross-varidation to determine lambda
        for (j in 1:length(alpha)){
          cv = cv.glmnet(X, Y, family = "binomial", alpha = alpha[j], nfolds=nfold)
          ind = match(cv$lambda.min, cv$lambda)
          cv.mat[j,2] = cv$lambda[ind]
          cv.mat[j,3] = cv$cvm[ind]
        }
        
        alpha.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),1]
        lambda.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),2]
        
        # Fit a GLM with elastic net regularization
        fit = glmnet(X.new, Y.new, family = "binomial", alpha = alpha.opt, lambda = lambda.opt)
        
        # Get model coefficients for glmnet
        coef = coef(fit)  
        
        if (iter == 1) {
          elast = as.matrix(coef)
        } else{
          elast = cbind(elast, as.matrix(coef))
        }
      } #end (iter)
      
      ### Get features for which coefficients are not 0
      elast = elast[-1, ] #get rid of intercept
      feature = rownames(elast)
      df = data.frame(elast, rowsum = rowSums(elast))
      ind.elast = which(df$rowsum !=0 ) #index of coefficients that are not zero 
      allFeatSel = as.character(feature[ind.elast])
      
      if (length(allFeatSel) != 0){
        ind = which(colnames(X) %in% allFeatSel)
        X.elast = as.matrix(X[,ind])
        colnames(X.elast) = allFeatSel
        # GLMs-NB Model
        pvalue.mat = get.glmNB(X.elast,Y) 
      } else{
        print("No significantly differentially abundant features!!!")
      }
      
      ### Summary of significantly differentially abundant features
      sig.mat = pvalue.mat[pvalue.mat[,"p.adj"] < threshold,] # Significantly differentially abundant features 
      sig = rownames(sig.mat)
      ind.sig = which(colnames(X) %in% sig)
      X.sig = t(as.matrix(X[,ind.sig]))
      X.sig.pop1 = X.sig[,1:npop1]
      X.sig.pop2 = X.sig[,(npop1+1):(npop1+npop2)]
      mean.pop1 = apply(X.sig.pop1, 1, mean)
      mean.pop2 = apply(X.sig.pop2, 1, mean)
      sd.pop1 = apply(X.sig.pop1, 1, sd)
      sd.pop2 = apply(X.sig.pop2, 1, sd)
      
      sig.df = data.frame(Annotation=as.character(rownames(X.sig)), mean.pop1,sd.pop1,mean.pop2, sd.pop2, 
                          sig.mat[,"p.val"], sig.mat[,"p.adj"], stringsAsFactors=F)
      names(sig.df)=c("Annotation", paste("mean_", names(table(Y))[1], sep=""), paste("sd_", names(table(Y))[1], sep=""), paste("mean_", names(table(Y))[2], sep=""), paste("sd_", names(table(Y))[2], sep=""), "p.val", "p.adj")                   
      
      rownames(sig.df) = NULL
      sigsort.df = sig.df[order(sig.df$p.adj),] 
      
    }
    
    alpha = seq(0.01,0.1,by=0.01)
    
    ### Data preprocessing and TMM normalization ###
    X = as.matrix(X)
    Y = factor(Y)
    X.tmm = TMMNorm(X,Y,TMM.option)  # TMM normalization (edgeR)
    X = t(X.tmm)
    
    ngroups = nlevels(Y)
    
    if (ngroups == 2) Twostage(X, Y) # Two group comparison

  }
  
  res <- tryCatch(TwoStage_Package(count_table,sample(outcome),TMM.option = 1),error = function(e) NULL)
  
  features <- data.frame(Feature = rownames(count_table))
  
  if(is.null(res)){
    output <- features
    output$pval <- 1
    output$pval.adj <- 1
  } else {
    colnames(res)[c(1,6:7)] <- c("Feature","pval","pval.adj")
    output <- merge(res,features, by="Feature", all = TRUE)
    output[is.na(output$pval),"pval"] <- 1
    output$pval.adj <- p.adjust(output$pval, method = p.adj)
  }
  output$Method <- "ENNB"
  return(output)
  
}



