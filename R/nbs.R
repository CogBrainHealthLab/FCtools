#' @title NBS
#'
#' @description Network-based statistics analysis
#'
#' @details This function implements the NBS analysis described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811910008852)}{Zalesky et al. (2010)}
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param nperm The number of permutations to generate the null distribution of network strengths. Set to 100 by default
#' @param nthread The number of CPU threads to use. Set to 1 by default
#' @param p the edge-wise threshold. Set to 0.001 by default
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
#' @export
############################################################################################################################
############################################################################################################################
NBS=function(model,contrast, FC_data, nperm=100, nthread=1, p=0.001)
{
  model=data.matrix(model)
  FC_data=data.matrix(FC_data)
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
  mod=.lm.fit(y = FC_data,x=data.matrix(cbind(1,model)))

  #define/init variables
  t.orig=extract.t(mod,colno+1)
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  tcrit=qt(p/2, NROW(model)-NCOL(model)-1, lower=FALSE)
  orig.clust=cluster.stat(t.orig,nnodes,tcrit)
  
  if(any(is.nan(orig.clust)){orig.clust=orig.clust[-which(is.nan(orig.clust[,1])),]} #if there are NaN values, which will happen for sparse SC matrices, they need to be recoded to 0s
  remove(mod)

  if(sum(orig.clust==0))  {cat(paste0("No significant networks are detected using the p<",p," edgewise threshold. You might want to use a more liberal threshold"))
  } else
  {
  ##permuted models
  #generating permutation sequences
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  for (perm in 1:nperm)  {permseq[,perm]=sample(1:NROW(model))}

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

  max.netstr=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("extract.t","cluster.stat"), .options.snow = opts)  %dopar%
    {
      #fitting permuted regression model and extracting max netstr in parallel streams
      mod.permuted=.lm.fit(y = FC_data,x=data.matrix(cbind(1,model))[permseq[,perm],])
      t.perm=extract.t(mod.permuted,colno+1)
      netstr=cluster.stat(t.perm,nnodes,tcrit)
      if(any(is.nan(netstr)){netstr=netstr[-which(is.nan(netstr[,1])),]} #if there are NaN values, which will happen for sparse SC matrices, they need to be recoded to 0s

      remove(t.perm,mod.permuted)

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
    returnobj=list(orig.clust,t.orig, tcrit,max.netstr)
    names(returnobj)=c("results","t.orig","tcrit","max.netstr")
    return(returnobj)
    suppressWarnings(closeAllConnections())
  }
}
############################################################################################################################
############################################################################################################################
#' @title extract.edges
#'
#' @description Generates edge-wise masks for calculating subject-level network strengths
#'
#' @details This function generates positive and negative masks (vectors of 1s and 0s), where 1s indicate a significant network-thresholded edge. These masks can then be used to perform a matrix multiplication with the vectorized FC matrices to object subject-level network strengths
#' @param NBS.obj A list object generated from an earlier `NBS()` analysis
#' @param network the network number (reported in the earlier NBS results) of the network to be masked. Set to 1 by default
#' @returns Returns a list object containing
#' \itemize{
#'  \item `clust.tstat` thresholded edge-wise t-statistics.Edges not belonging to this cluster will be zeroed.
#'  \item `pos.edges` A vector of 1s and 0s indicating the significant network-thresholded positive edges.
#'  \item `neg.edges` A vector of -1s and 0s indicating the significant network-thresholded negative edges.
#'  \item `pos.mask` A vector of 1s and 0s indicating the significant network-thresholded positive edges.
#'  \item `neg.mask` A vector of 1s and 0s indicating the significant network-thresholded negative edges.
#'}
#' @examples
#' \dontrun{
#' extract.edges(model1,network=1)
#' }
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export

extract.edges=function(NBS.obj,network=1)
{
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -length(NBS.obj$t.orig))) / (2 * 0.5)
  ##recode all p="<0.**" into 0 for subsequent thresholding
  if(is.character(NBS.obj$results[,4]))
  {
    NBS.obj$results[,4]=suppressWarnings(as.numeric(NBS.obj$results[,4]))
    NBS.obj$results[,4][is.na(NBS.obj$results[,4])]=0
  }
  if(is.character(NBS.obj$results[,5]))
  {
    NBS.obj$results[,5]=suppressWarnings(as.numeric(NBS.obj$results[,5]))
    NBS.obj$results[,5][is.na(NBS.obj$results[,5])]=0
  }

  ##thresholding tstats
  tstat.thresholded=NBS.obj$t.orig
  tstat.thresholded[abs(NBS.obj$t.orig)<NBS.obj$tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1

  ##reshaping 1D tstat vector to 2D matrices
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)

  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded-(tstat.thresholded.bin*NBS.obj$tcrit) ## subtracting tcrit values to be consist with NBR::nbr_lm()
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin

  ##clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  idx=which(com$membership==which(com$csize>2)[network])

  ##masking out edges from other networks
  FC_mat.mask=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.mask[idx,idx]=1
  mask=FC_mat.mask[upper.tri(FC_mat.mask,diag = F)]
  clust.tstat=tstat.thresholded*mask

    #positive mask
    clust.pos.mask=clust.tstat
    clust.pos.mask[clust.pos.mask>0]=1
    clust.pos.mask[clust.pos.mask<0]=0
  
    #negative mask
    clust.neg.mask=clust.tstat
    clust.neg.mask[clust.neg.mask>0]=0
    clust.neg.mask[clust.neg.mask<0]=-1

  ##objects to return
  returnobj=list(as.numeric(clust.tstat),as.numeric(clust.pos.mask),as.numeric(clust.neg.mask),as.numeric(clust.pos.mask),abs(as.numeric(clust.neg.mask)))
  names(returnobj)=c("clust.tstat","pos.edges","neg.edges","pos.mask","neg.mask")
  return(returnobj)
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

### extract t-stats efficiently
extract.t=function(mod,row)
{
  p = mod$rank
  df.residual=NROW(mod$residuals)-NROW(mod$coefficients)
  rdf = df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr[p1, p1, drop = FALSE])
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se
  return(tval)
}
