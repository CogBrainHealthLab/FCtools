#' @title network_lm
#' @description mass univariate linear regression at the network level
#'
#' @details This function first summarizes the FC edges into their respective networks and then carry out mass univariate linear regression analyses on each of the network to network connection
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param threshold.method method for correcting for multiple tests. set to `fdr` by default
#' @param perm If set to `TRUE`, p values will be calculated using a permutation approach by shuffling subjects' labels, before correcting for FDR. Set to `TRUE` by default
#' @param nperm number of permutations to use if `perm=TRUE`.
#' @returns Returns a data.frame object with `coef` and corrected `p` values
#'
#' @examples
#' \dontrun{
#' model1=network_lm(model,contrast, FC_data)
#' }
#' @export
############################################################################################################################
############################################################################################################################
network_lm=function(model,contrast, FC_data, threshold.method="fdr",perm=TRUE, nperm=1000)
{
  ##checks
  #check if nrow is consistent for model and FC_data
  if(NROW(FC_data)!=NROW(model))  {stop(paste("The number of rows for FC_data (",NROW(FC_data),") and model (",NROW(model),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(model)==FALSE)
  if(length(idxF)>0)
  {
    cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    model=model[-idxF,]
    contrast=contrast[-idxF]
    FC_data=FC_data[-idxF,]
  }
  
  #check contrast
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {stop("contrast is not contained within model")}
      
      if(class(contrast) != "integer" & class(contrast) != "numeric")
      {
        if(identical(contrast,model[,colno]))  {break}
      } else
      {
        if(identical(as.numeric(contrast),as.numeric(model[,colno])))  {break}
      }
    }
  }  else
  {
    if(class(contrast) != "integer" & class(contrast) != "numeric")
    {
      if(identical(contrast,model))  {colno=1}
      else  {stop("contrast is not contained within model")}
    } else
    {
      if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
      else  {stop("contrast is not contained within model")}
    }
  }
  
  #check categorical variable
  if(NCOL(model)>1)
  {
    for (column in 1:NCOL(model))
    {
      if(class(model[,column]) != "integer" & class(model[,column]) != "numeric")
      {
        if(length(unique(model[,column]))==2)
        {
          cat(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model[,column]==unique(model[,column])[2]]=1
          model[,column]=recode
          contrast=model[,colno]
        } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }
    }
  } else
  {
    if (!suppressWarnings(all(!is.na(as.numeric(as.character(model))))))
    {
      if(length(unique(model))==2)
      {
        cat(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(model))
        recode[model==unique(model)[2]]=1
        model=recode
        contrast=model
      } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }
  }
  
  
  
  model=data.matrix(model)
  FC_data=data.matrix(FC_data)
  
  #collinearity check
  if(NCOL(model)>1)
  {
    cormat=cor(model,use = "pairwise.complete.obs")
    cormat.0=cormat
    cormat.0[cormat.0==1]=NA
    if(max(abs(cormat.0),na.rm = TRUE) >0.7)
    {
      warning(paste("correlations among variables in model are observed to be as high as ",round(max(abs(cormat.0),na.rm = TRUE),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...",sep=""))
    }
  }
  #fit model
  FNC_data=scale(FCtools::edges_to_networks(FC_data))
  model.std=data.matrix(scale(model))
  mod.fitted=.lm.fit(y = FNC_data,x=model.std)
  
  #compute pvalues
  
  n=nrow(FNC_data)
  p=ncol(model.std)
  Nedges=ncol(FNC_data)
  
  XtX_inv <- solve(t(model.std) %*% model.std)
  results.lm <- lapply(1:ncol(FNC_data), function(j) {
    resid_j <- mod.fitted$residuals[, j]
    sigma2_j <- sum(resid_j^2) / (n - p)
    se_j <- sqrt(diag(XtX_inv) * sigma2_j)
    t_vals_j <- mod.fitted$coefficients[, j] / se_j
    p_vals_j <- 2 * pt(abs(t_vals_j), df = n - p, lower.tail = FALSE)})
  
  if(perm==TRUE)
  { set.seed(123)
    #permutation
    coef.unperm=mod.fitted$coefficients[colno,]
    coef.perm=matrix(NA,nrow=nperm, ncol=Nedges)
    for(iter in 1:nperm)
    {
      mod.fitted=.lm.fit(y = FNC_data[sample(n),],x=model.std)
      coef.perm[iter,]=as.numeric(mod.fitted$coefficients[colno,])
    }
    
    p.perm=rep(NA,Nedges)
    #calculate permutation p values
    for(connection in 1:Nedges)
    {
      if(coef.unperm[connection]<0)
      {
        if(length(which(coef.perm[,connection]<coef.unperm[connection]))>0)
        {
          p.perm[connection]=length(which(coef.perm[,connection]<coef.unperm[connection]))/nperm  
        } else
        {
          p.perm[connection]=1/nperm
        }  
      } else
      {
        if(length(which(coef.perm[,connection]>coef.unperm[connection]))>0)
        {
          p.perm[connection]=length(which(coef.perm[,connection]>coef.unperm[connection]))/nperm  
        } else
        {
          p.perm[connection]=1/nperm
        }  
      }
    }
    results=data.frame(coef=coef.unperm,
                       p.thresholded=p.adjust(p.perm,method = threshold.method))
  } else
  {
    results=data.frame(coef=coef.unperm,
                       p.thresholded=p.adjust(t(matrix(unlist(results.lm), nrow=ncol(model)))[,colno],method = threshold.method))
  }
  rownames(results)=colnames(FNC_data)
  return(results)
}