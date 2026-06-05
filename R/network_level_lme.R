#' @title network_lme
#' @description mass univariate linear mixed effects analysis at the network level
#'
#' @details This function first summarizes the FC edges into their respective networks and then carry out mass univariate linear mixed effect analyses on each of the network to network connection
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param threshold.method method for correcting for multiple tests. set to `fdr` by default
#' @returns Returns a data.frame object with `coef` and corrected `p` values
#'
#' @examples
#' \dontrun{
#' model1=network_lme(model,contrast, FC_data)
#' }
#' @export
############################################################################################################################
############################################################################################################################
network_lme=function(model,contrast,random, FC_data,threshold.method="fdr")
{
  ##checks
  #check random variable
  if(missing(random))   {stop("The 'random' parameter has to be specified")}
  else  {random=match(random,unique(random))}
  
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
    random=random[-idxF]
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
        if(identical(contrast,model[,colno]))  {break}
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
      if(identical(contrast,model))  {colno=1}
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
  
  
  
  model=scale(data.matrix(model))
  # formula_dataset.copy=formula_dataset
  FNC_data=scale(edges_to_networks(data.matrix(FC_data)))
  Nedges=NCOL(FNC_data)
  
  coef=rep(NA,Nedges)
  p=rep(NA,Nedges)
  
  for(connection in 1:Nedges)
  {
    results.connection=lmefast.p(FNC_data[,connection], model, random, contrast=colno+1)
    coef[connection]=results.connection[1]
    p[connection]=results.connection[2]
  }
  results=data.frame(coef=coef,p.thresholded=p.adjust(p,method=threshold.method))
  rownames(results)=colnames(FNC_data)
    
  return(results)
}


lmefast.p=function(Y,x, id, tol = 1e-07, ranef = TRUE, maxiters = 200, contrast)
{
  # Force numeric matrices
  x <- data.matrix(x)
  id <- as.integer(as.factor(id))
  Y <- as.matrix(Y)
  k <- ncol(Y)
  
  # Apply rint.reg efficiently across Y columns
  beta.p<- sapply(1:k, function(j) {
    mod <- Rfast::rint.reg(Y[, j], x, id, tol = tol, ranef = ranef, maxiters = maxiters)
    c(mod$be[contrast],2 * (1 - pnorm(abs(mod$be[contrast] / mod$se[contrast]))))
  })
  
  return(beta.p)
}