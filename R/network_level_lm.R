#' @title network_lm
#' @description mass univariate linear regression at the network level
#'
#' @details This function implements the NBS analysis described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811910008852)}{Zalesky et al. (2010)}
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param threshold.method method for correcting for multiple tests. set to `fdr` by default
#' @returns Returns a data.frame object with `coef` and corrected `p` values
#'
#' @examples
#' \dontrun{
#' model1=network_lm(model,contrast, FC_data)
#' }
#' @export
############################################################################################################################
############################################################################################################################
network_lm=function(model,contrast, FC_data, threshold.method="fdr")
{
  ##checks

  #check if nrow is consistent for model and FC_data
  if(NROW(FC_data)!=NROW(model))  {stop(paste("The number of rows for FC_data (",NROW(FC_data),") and model (",NROW(model),") are not the same",sep=""))}

  #incomplete data check
  idxF=which(complete.cases(model)==F)
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
    if(max(abs(cormat.0),na.rm = T) >0.5)
    {
      warning(paste("correlations among variables in model are observed to be as high as ",round(max(abs(cormat.0),na.rm = T),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...",sep=""))
    }
  }
  #fit model
  FNC_data=scale(edges_to_networks(FC_data))
  model.std=data.matrix(scale(model))
  mod.fitted=.lm.fit(y = FNC_data,x=model.std)

  #compute pvalues
    
  n=nrow(FNC_data)
  p=ncol(model.std)
  
  XtX_inv <- solve(t(model.std) %*% model.std)
  results.lm <- lapply(1:ncol(FNC_data), function(j) {
  resid_j <- mod.fitted$residuals[, j]
  sigma2_j <- sum(resid_j^2) / (n - p)
  se_j <- sqrt(diag(XtX_inv) * sigma2_j)
  t_vals_j <- mod.fitted$coefficients[, j] / se_j
  p_vals_j <- 2 * pt(abs(t_vals_j), df = n - p, lower.tail = FALSE)})
    
  results=data.frame(coef=t(mod.fitted$coefficients)[,colno],p.thresholded=p.adjust(t(matrix(unlist(results.lm), nrow=2))[,colno],method = threshold.method))
  rownames(results)=colnames(FNC_data)
    
  return(results)
}

