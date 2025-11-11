#' @title NBS
#'
#' @description Network-based statistics analysis
#'
#' @details This function implements the NBS analysis described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811910008852)}{Zalesky et al. (2010)}
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param nperm The number of permutations to generate the null distribution of network strengths. Set to 100 by default
#' @param nthread The number of CPU threads to use. Set to 1 by default
#' @param p the edge-wise threshold. Set to 0.001 by default
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#' @returns Returns a list object containing
#' \itemize{
#'  \item `results` Edge- and network-wise results in a data.frame object
#'  \item `t.orig` Edge-wise t-stats
#'  \item `tcrit` The critical t-value
#'  \item `max.netstr` A vector containing the null distribution of the permuted network strengths
#'}
#' @examples
#' \dontrun{
#' model1=NBS(model,contrast, FC_data, nperm=1000, nthread=8, p=0.001)
#' }
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom Rfast rint.reg
#' @export
############################################################################################################################
############################################################################################################################
NBS_lme=function(model,contrast,random, FC_data, nperm=100, nthread=1, p=0.001,perm_type="row")
{
  model=data.matrix(model)
  FC_data=data.matrix(FC_data)
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

  ##unpermuted model
  z.orig=lmefast(FC_data, model, random,contrast=colno+1)
  tcrit=qnorm(p/2, mean = 0, sd = 1, lower.tail = FALSE)
  
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  orig.clust=cluster.stat(z.orig,nnodes, tcrit)
  
  if(any(is.nan(orig.clust))){orig.clust=orig.clust[-which(is.nan(orig.clust[,1])),]} #if there are NaN values, which will happen for sparse SC matrices, they need to be recoded to 0s

  if(sum(orig.clust==0))  {cat(paste0("No significant networks are detected using the p<",p," edgewise threshold. You might want to use a more liberal threshold"))
  } else
  {
  ##permuted models
  #generating permutation sequences
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
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


  start=Sys.time()
  cat("\nEstimating permuted network strengths...\n")

  max.netstr=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("cluster.stat","lmefast"),.packages = "Rfast", .options.snow = opts)  %dopar%
    {
      #fitting permuted regression model and extracting max netstr in parallel streams
      z.perm=lmefast(FC_data, model[permseq[,perm],], random,contrast=colno+1)
      netstr=cluster.stat(z.perm,nnodes,tcrit)
      if(any(is.nan(netstr))) {netstr=netstr[-which(is.nan(netstr[,1])),]} #if there are NaN values, which will happen for sparse SC matrices, they need to be recoded to 0s

      remove(z.perm)

      if(length(netstr)>2)  {max.netstr=c(max(netstr[,1]),max(netstr[,2]))}
      else {max.netstr=netstr} #if only one row of results is obtained, there is no need to use the max() function

      return(max.netstr)
    }
  end=Sys.time()
  cat(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))

  ##processing results
    #saving cluster-related results into a data.frame object
    orig.clust=data.frame(orig.clust)
    orig.clust$p.unweighted=NA
    orig.clust$p.weighted=NA

    #thresholding clusters using permuted null distribution
    for(row in 1:nrow(orig.clust))
    {
      orig.clust[row,3]=sum(max.netstr[,1] > orig.clust[row,1])/nperm
      orig.clust[row,4]=sum(max.netstr[,2] > orig.clust[row,2])/nperm
    }

    #formatting results table
    orig.clust[,c(3,4)][orig.clust[,c(3,4)]==0]=paste("<",1/nperm,sep="") #if p=0
    orig.clust=cbind(c(1:nrow(orig.clust)),orig.clust)
    colnames(orig.clust)=c("network","strength.unweighted","strength.weighted","p.unweighted","p.weighted")
  
    #objects to return
    returnobj=list(orig.clust,z.orig, tcrit,max.netstr)
    names(returnobj)=c("results","t.orig","tcrit","max.netstr")
    return(returnobj)
    suppressWarnings(closeAllConnections())
  }
}

############################################################################################################################
############################################################################################################################
### extract cluster-stats from graphs
cluster.stat=function(data,nnodes,tcrit)
{
  ##thresholding
  tstat.thresholded=data
  tstat.thresholded[abs(data)<tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1
  
  ##setting up FCmatrices
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)
  
  ##thresholding
  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin
  FC_mat.weighted=abs(FC_mat.weighted)-(FC_mat.unweighted*tcrit)
  
  #clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  
  #count edges in clusters
  if(length(which(com$csize>2)>0)) #proceed only if there is at least one cluster with 3 nodes (i.e 2 edges). Isolated/unconnected edges are removed
  {
    cluster.idx=which(com$csize>2)
    clust.results=matrix(NA,nrow=length(cluster.idx), ncol=2)
    
    for (cluster.no in 1:length(cluster.idx))
    {
      idx=which(com$membership==cluster.idx[cluster.no])
      clust.results[cluster.no,1]=strength.unweighted=sum(FC_mat.unweighted[idx,idx])
      clust.results[cluster.no,2]=strength.weighted=sum(FC_mat.weighted[idx,idx])
    }
  } else  { clust.results=c(0,0)} #if no clusters (with 3 nodes) are detected
  return(clust.results)
}
############################################################################################################################
############################################################################################################################

## permutation functions for random subject effects
## Paired/grouped data points are first shuffled within subjects, then these pairs/groups are shuffled between subjects
perm_within_between=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count))>0)
    {
      sub.id=as.numeric(which(table(random)==count))
      if(length(sub.id)>1)
      {
        ##between group shuffling
        recode.vec=sample(sub.id)
        vec.idx=1
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=sample(which(random==recode.vec[vec.idx])) ##sample— within subject shuffling
          vec.idx=vec.idx+1
        }   
        remove(vec.idx,recode.vec)  
      } else 
      {
        ##if only one subject has a certain count, between subject shuffling will not be possible, only within-subject shuffling will be carried out
        perm.idx[which(random==sub.id)]=sample(which(random==sub.id)) ##sample— within subject shuffling
      }
    }
  }
  ##for subjects with a single measurement
  sub.idx=which(is.na(perm.idx))
  if(length(sub.idx)>1)
  {
    perm.idx[sub.idx]=sample(sub.idx)  
  } else 
  {
    perm.idx[sub.idx]=sub.idx
  }
  return(perm.idx)
}

## Paired/grouped data points are shuffled within subjects, order of subjects in the dataset remains unchanged
perm_within=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count)>0))
    {
      sub.id=as.numeric(which(table(random)==count))
      for(sub in sub.id)
      {
        perm.idx[which(random==sub)]=sample(which(random==sub))
      }  
    }
  }
  return(perm.idx)
}

## Paired/grouped data points are shuffled between subjects, order of data points within subjects remains unchanged.
perm_between=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count))>0)
    {
      sub.id=as.numeric(which(table(random)==count))
      if(length(sub.id)>1)
      {
        ##between group shuffling
        recode.vec=sample(sub.id)
        vec.idx=1
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=which(random==recode.vec[vec.idx])
          vec.idx=vec.idx+1
        }   
        remove(vec.idx,recode.vec)  
      }
    }
  }
  ##for subjects with a single measurement
  sub.idx=which(is.na(perm.idx))
  if(length(sub.idx)>1)
  {
    perm.idx[sub.idx]=sample(sub.idx)  
  } else 
  {
    perm.idx[sub.idx]=sub.idx
  }
  return(perm.idx)
}

############################################################################################################################
############################################################################################################################

lmefast=function(Y,x, id, tol = 1e-07, ranef = TRUE, maxiters = 200, contrast)
{
  # Force numeric matrices
  x <- data.matrix(x)
  id <- as.integer(as.factor(id))
  Y <- as.matrix(Y)
  k <- ncol(Y)
  
  # Apply rint.reg efficiently across Y columns
  tstats <- sapply(1:k, function(j) {
    mod <- Rfast::rint.reg(Y[, j], x, id, tol = tol, ranef = ranef, maxiters = maxiters)
    mod$be[contrast] / mod$se[contrast]
  })
  
  return(tstats)
}

