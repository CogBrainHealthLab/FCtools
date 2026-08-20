#' @title network_lme
#' @description mass univariate linear mixed effects analysis at the network level
#'
#' @details This function first summarizes the FC edges into their respective networks and then carry out mass univariate linear mixed effect analyses on each of the network to network connection
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param perm If set to `TRUE`, p values will be calculated using a permutation approach by shuffling subjects' labels, before correcting for FDR. Set to `TRUE` by default
#' @param nperm number of permutations to use if `perm=TRUE`.
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#' @param threshold.method method for correcting for multiple tests. set to `fdr` by default
#' @returns Returns a data.frame object with `coef` and corrected `p` values
#'
#' @export
############################################################################################################################
############################################################################################################################
network_lme=function(model,contrast,random, FC_data,threshold.method="fdr",perm=TRUE, nperm=1000,perm_type="within_between",nthread=4)
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
    random=random[-idxF]
  }
  
  #check random variable
  if(missing(random))   {stop("The 'random' parameter has to be specified")}
  else  {random=match(random,unique(random))}
  
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
  n=NROW(FNC_data)
  coef.unperm=rep(NA,Nedges)
  p=rep(NA,Nedges)
  
  for(connection in 1:Nedges)
  {
    results.connection=lmefast.p(FNC_data[,connection], model, random, contrast=colno+1)
    coef.unperm[connection]=results.connection[1]
    p[connection]=results.connection[2]
  }
  
  if(perm==TRUE)
  {
    start=Sys.time()
    #generating permutation sequences
    set.seed(123)
    permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
    
    if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(as.numeric(random))}} 
    else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
    else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
    else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
    
    #activate parallel processing
    unregister_dopar = function() {
      env = foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    
    cl=parallel::makeCluster(nthread)
    doParallel::registerDoParallel(nthread)
    `%dopar%` = foreach::`%dopar%`
    
    #progress bar
    doSNOW::registerDoSNOW(cl)
    pb=txtProgressBar(max = nperm, style = 3)
    progress=function(n) setTxtProgressBar(pb, n)
    opts=list(progress = progress)
    
    coef.perm=matrix(NA,nrow=nperm, ncol=Nedges)
    pb = txtProgressBar(min = 0, max = nperm, style = 3) 
    cat("\nEstimating permuted network strengths...\n")
    coef.perm.all=foreach::foreach(iter=1:nperm, .combine="rbind",.export="lmefast.p",.packages = "Rfast", .options.snow = opts)  %dopar%
      {
        coef.perm=matrix(NA,nrow=1, ncol=Nedges)
        for(connection in 1:Nedges)
        {
          results.connection=lmefast.p(FNC_data[permseq[,iter],connection], model, random, contrast=colno+1,p=FALSE)
          coef.perm[1,connection]=results.connection[1]
        }
        return(coef.perm)
      }
    closeAllConnections()
    p.perm=rep(NA,Nedges)
    #calculate permutation p values
    for(connection in 1:Nedges)
    {
      if(coef.unperm[connection]<0)
      {
        if(length(which(coef.perm.all[,connection]<coef.unperm[connection]))>0)
        {
          p.perm[connection]=length(which(coef.perm.all[,connection]<coef.unperm[connection]))/nperm  
        } else
        {
          p.perm[connection]=1/nperm
        }  
      } else
      {
        if(length(which(coef.perm.all[,connection]>coef.unperm[connection]))>0)
        {
          p.perm[connection]=length(which(coef.perm.all[,connection]>coef.unperm[connection]))/nperm  
        } else
        {
          p.perm[connection]=1/nperm
        }  
      }
    }
    results=data.frame(coef=coef.unperm,
                       p.thresholded=p.adjust(p.perm,method = threshold.method))
    end=Sys.time()
    cat(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  } else
  {
    results=data.frame(coef=coef.unperm,p.thresholded=p.adjust(p,method=threshold.method))
  }

  rownames(results)=colnames(FNC_data)
  return(results)
}


lmefast.p=function(Y,x, id, tol = 1e-07, ranef = TRUE, maxiters = 200, contrast,p=TRUE)
{
  # Force numeric matrices
  x <- data.matrix(x)
  id <- as.integer(as.factor(id))
  Y <- as.matrix(Y)
  k <- ncol(Y)
  
  # Apply rint.reg efficiently across Y columns
  if(p==TRUE)
  {
    beta.p<- sapply(1:k, function(j) {
      mod <- Rfast::rint.reg(Y[, j], x, id, tol = tol, ranef = ranef, maxiters = maxiters)
      c(mod$be[contrast],2 * (1 - pnorm(abs(mod$be[contrast] / mod$se[contrast]))))
    })  
  }
  
  if(p==FALSE)
  {
    beta.p<- sapply(1:k, function(j) {
      mod <- Rfast::rint.reg(Y[, j], x, id, tol = tol, ranef = ranef, maxiters = maxiters)
      mod$be[contrast]
    })
    
  }
  return(beta.p)
}
