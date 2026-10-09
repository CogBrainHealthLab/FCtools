#' @title NBS
#'
#' @description Network-based statistics analysis
#'
#' @details This function implements the NBS analysis described in
#' Zalesky et al. (2010) \doi{10.1016/j.neuroimage.2010.06.041}
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param nperm The number of permutations to generate the null distribution of network strengths. Set to 100 by default
#' @param nthread The number of CPU threads to use. Set to 1 by default
#' @param p the edge-wise threshold. Set to 0.001 by default
#' @returns A list object containing
#' \itemize{
#'  \item `results` Edge- and network-wise results in a data.frame object
#'  \item `t.orig` Edge-wise t-stats
#'  \item `tcrit` The critical t-value
#'  \item `max.netstr` A vector containing the null distribution of the permuted network strengths
#'}
#' @examples
#' demomat=get('demomat')[,1:7021] 
#' contrast=c(1,1,2,2)
#' random=c('sub1','sub2','sub3','sub4')
#' model1=NBS(model=contrast, 
#'            contrast=contrast, 
#'            FC_data=demomat, 
#'            nperm=2, 
#'            nthread=2, 
#'            p=0.001)
#' @useDynLib FCtools, .registration = TRUE
#' @importFrom mori share
#' @importFrom Rcpp evalCpp
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils getFromNamespace txtProgressBar setTxtProgressBar
#' @importFrom stats complete.cases cor qt
#' @export
############################################################################################################################
############################################################################################################################
NBS=function(model,contrast, FC_data, nperm=100, nthread=1, p=0.001)
{
  ##checks

  #check if nrow is consistent for model and FC_data
  if(NROW(FC_data)!=NROW(model))  {stop(paste("The number of rows for FC_data (",NROW(FC_data),") and model (",NROW(model),") are not the same",sep=""))}

  #incomplete data check
  idxF=which(stats::complete.cases(model)==FALSE)
  if(length(idxF)>0)
  {
    warning(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    if(!is.null(ncol(model))) {model=model[-idxF,]} else {model=model[-idxF]}
    contrast=contrast[-idxF]
    FC_data=FC_data[-idxF,]
  }

  #check contrast
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {stop("contrast is not contained within model")}

      if(!inherits(contrast,'integer') & !inherits(contrast,'numeric'))
      {
        if(identical(contrast,model[,colno]))  {break}
      } else
      {
        if(identical(suppressWarnings(as.numeric(contrast)),
                     suppressWarnings(as.numeric(model[,colno]))))  {break}
      }
    }
  }  else
  {
    if(!inherits(contrast,'integer') & !inherits(contrast,'numeric'))
    {
      if(identical(contrast,model))  {colno=1}
      else  {stop("contrast is not contained within model")}
    } else
    {
      if(identical(suppressWarnings(as.numeric(contrast)),
                   suppressWarnings(as.numeric(model))))  {colno=1}
      else  {stop("contrast is not contained within model")}
    }
  }

  #check categorical variable
  if(NCOL(model)>1)
  {
    for (column in 1:NCOL(model))
    {
      if(!inherits(model[,column],'integer') & !inherits(model[,column],'numeric'))
      {
        if(length(unique(model[,column]))==2)
        {
          message(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))

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
        message(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))

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
    cormat=stats::cor(model,use = "pairwise.complete.obs")
    cormat.0=cormat
    cormat.0[cormat.0==1]=NA
    if(max(abs(cormat.0),na.rm = TRUE) >0.5)
    {
      warning(paste("correlations among variables in model are observed to be as high as ",round(max(abs(cormat.0),na.rm = TRUE),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...",sep=""))
    }
  }

  ## Unpermuted model
  X <- data.matrix(cbind(Intercept = 1, model))
  tested <- colno + 1L
  
  # Constant edges, including structural zeros, cannot supply a t-statistic.
  constant.edges <- vapply(
    seq_len(ncol(FC_data)),
    function(j) all(FC_data[, j] == FC_data[1L, j]),
    logical(1)
  )
  
  mod <- lm_fast(X, FC_data)
  t.orig <- as.numeric(mod[tested, ])
  t.orig[constant.edges] <- 0

  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  tcrit=stats::qt(p/2, NROW(model)-NCOL(model)-1, lower=FALSE)
  orig.clust=cluster.stat(t.orig,nnodes,tcrit)
  
  if(any(is.nan(orig.clust))){orig.clust=orig.clust[-which(is.nan(orig.clust[,1])),]} #if there are NaN values, which will happen for sparse SC matrices, they need to be recoded to 0s
  remove(mod)
  suppressWarnings(closeAllConnections())
  if(sum(orig.clust==0))  {message(paste0("No significant networks are detected using the p<",p," edgewise threshold. You might want to use a more liberal threshold"))
    
  } else
  {

    # Freedman-Lane reduced model:
    # retain the intercept and nuisance predictors; remove the tested predictor.
    Z <- X[, -tested, drop = FALSE]
    
    reduced <- lm_fast(Z, FC_data, return_coefficients = TRUE)
    fitted0 <- Z %*% reduced$coefficients
    residuals0 <- FC_data - fitted0
    
    # Generate row permutations in the main process.
    # Each permutation uses the same row order across all edges.
    permseq <- matrix(0L, nrow = nrow(X), ncol = nperm)
    
    for (perm in seq_len(nperm)) {
      permseq[, perm] <- sample.int(nrow(X))
    }
    
    # Shared read-only inputs for workers on the same machine.
    shared <- mori::share(list(
      X = X,
      fitted0 = fitted0,
      residuals0 = residuals0,
      permseq = permseq
    ))
    
    # Keep 'shared' alive until the workers have finished.
    rm(Z, reduced, fitted0, residuals0, permseq, X)
    
  ##permuted models
    # Local PSOCK workers work on Windows, macOS and Linux.
    cl <- parallel::makeCluster(nthread, type = "PSOCK")
    pb <- NULL
    
    on.exit({
      if (!is.null(pb)) try(close(pb), silent = TRUE)
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    # Load mori before workers deserialize the shared objects.
    parallel::clusterEvalQ(cl, {
      loadNamespace("mori")
      loadNamespace("FCtools")
      NULL
    })
    
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    
    pb <- utils::txtProgressBar(max = nperm, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  
  #safe connection clean up at the end of the run
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  start=Sys.time()
  message("\nEstimating permuted network strengths...\n")

  max.netstr <- foreach::foreach(
    perm = seq_len(nperm),
    .combine = "rbind",
    .packages = c("FCtools", "mori"),
    .export = c("cluster.stat"),
    .options.snow = opts
  ) %dopar% {
    
    # Load the installed wrapper and its package-native C++ binding.
    # Freedman-Lane:
    # Y* = nuisance fitted values + permuted nuisance-model residuals.
    Yperm <- shared$fitted0 +
      shared$residuals0[
        shared$permseq[, perm],
        ,
        drop = FALSE
      ]
    
    # The full design matrix remains unchanged.
    t.perm <- as.numeric(
      lm_fast(shared$X, Yperm)[tested, ]
    )
    
    t.perm[constant.edges] <- 0
    
    netstr <- cluster.stat(t.perm, nnodes, tcrit)
    
    if (is.null(dim(netstr))) {
      c(0, 0)  # No qualifying network
    } else {
      c(max(netstr[, 1L]),max(netstr[, 2L])
      )
    }
  }
  
  # Preserve two-dimensional indexing even when nperm = 1.
  max.netstr <- matrix(max.netstr, ncol = 2L)
  end=Sys.time()
  message(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))

  ##processing results
    #saving cluster-related results into a data.frame object
    orig.clust=data.frame(orig.clust)
    orig.clust$p.unweighted=NA
    orig.clust$p.weighted=NA

    #thresholding clusters using permuted null distribution
    for(row in 1:nrow(orig.clust))
    {
      orig.clust[row, 3] <-round((1 + sum(max.netstr[, 1] >= orig.clust[row, 1])) / (nperm + 1),4)
      orig.clust[row, 4] <-round((1 + sum(max.netstr[, 2] >= orig.clust[row, 2])) / (nperm + 1),4)
    }

    #formatting results table
    orig.clust=cbind(c(1:nrow(orig.clust)),orig.clust)
    colnames(orig.clust)=c("network","strength.unweighted","strength.weighted","p.unweighted","p.weighted")
  
    #objects to return
    returnobj=list(orig.clust,t.orig, tcrit,max.netstr)
    names(returnobj)=c("results","t.orig","tcrit","max.netstr")
    return(returnobj)
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
#' demomat=get('demomat')
#' contrast=c(1,1,2,2)
#' random=c('sub1','sub2','sub3','sub4')
#' model1=NBS(model=contrast, contrast=contrast, FC_data=demomat, nperm=2, nthread=1, p=0.001)
#' 
#' edges=extract.edges(model1,network=1)
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

  FC_mat.weighted[upper.tri(FC_mat.weighted,diag =FALSE)]=tstat.thresholded-(tstat.thresholded.bin*NBS.obj$tcrit) ## subtracting tcrit values to be consist with NBR::nbr_lm()
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag =FALSE)]=tstat.thresholded.bin

  ##clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='max', weighted=NULL))
  idx=which(com$membership==which(com$csize>2)[network])

  ##masking out edges from other networks
  FC_mat.mask=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.mask[idx,idx]=1
  mask=FC_mat.mask[upper.tri(FC_mat.mask,diag =FALSE)]
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
### Extract cluster statistics using C++
cluster.stat <- function(data, nnodes, tcrit) {
  cluster_stat_cpp(
    data = as.numeric(data),
    nnodes = as.double(nnodes),
    tcrit = as.double(tcrit)
  )
}

