#' @title NBS_lme
#'
#' @description Network-based statistics analysis
#'
#' @details This function implements the NBS analysis described in
#' Zalesky et al. (2010) \doi{10.1016/j.neuroimage.2010.06.041}, using linear mixed effects models.
#' @param model A data.frame or matrix containing all the predictors in the model
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param contrast The predictor of interest. The edge- and network-wise statistics will only be estimated for this predictor
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param nperm The number of permutations to generate the null distribution of network strengths. Set to 100 by default
#' @param nthread The number of CPU threads to use. Set to 1 by default
#' @param p the edge-wise threshold. Set to 0.001 by default
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#' @returns A list object containing
#' \itemize{
#'  \item `results` Edge- and network-wise results in a data.frame object
#'  \item `t.orig` Edge-wise t-stats
#'  \item `tcrit` The critical t-value
#'  \item `max.netstr` A vector containing the null distribution of the permuted network strengths
#'}
#' @examples
#' demomat=get('demomat')[,1:4005] 
#' contrast=c(1,1,2,2)
#' random=c('sub1','sub2','sub3','sub4')
#' model1=NBS_lme(model=contrast, 
#'                contrast=contrast, 
#'                random=random, 
#'                FC_data=demomat, 
#'                nperm=1,
#'                nthread=1, 
#'                p=0.001)
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
NBS_lme=function(model,contrast,random, FC_data, nperm=100, nthread=1, p=0.001,perm_type="row")
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
      
      if(!inherits(contrast, "integer") & !inherits(contrast, "numeric"))
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
    if(!inherits(contrast, "integer") & !inherits(contrast, "numeric"))
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

  ## Unpermuted random-intercept model
  X <- data.matrix(cbind(Intercept = 1, model))
  tested <- colno + 1L
  nedge <- ncol(FC_data)
  
  # Exclude constant edges from variance-component estimation.
  # Restore these edges as zero statistics in the full edge vector.
  active <- which(vapply(
    seq_len(nedge),
    function(j) any(FC_data[, j] != FC_data[1L, j]),
    logical(1)
  ))
  
  Y <- FC_data[, active, drop = FALSE]
  fit.orig <- lme_fast(Y = Y,X = X,id = random,add_intercept = FALSE)
  z.orig <- numeric(nedge)
  z.orig[active] <- as.numeric(fit.orig$t_stat[tested, ])

  # Retain the original normal-reference cluster-forming threshold.
  # Despite the variable name, z.orig contains coefficient/SE statistics.
  tcrit <- stats::qnorm(p / 2, lower.tail = FALSE)
  nnodes <- (1 + sqrt(1 + 8 * nedge)) / 2
  orig.clust <- cluster.stat(z.orig, nnodes, tcrit)
  
  rm(fit.orig)

  if (is.null(dim(orig.clust))) {
    message(
      "No qualifying networks detected at the p < ", p,
      " edgewise threshold."
    )
    
    return(list(
      results = data.frame(
        network = integer(),
        strength.unweighted = numeric(),
        strength.weighted = numeric(),
        p.unweighted = numeric(),
        p.weighted = numeric()
      ),
      t.orig = z.orig,
      tcrit = tcrit,
      max.netstr = matrix(numeric(), nrow = 0L, ncol = 2L)
    ))
  } else
  {
  ##permuted models
  #generating permutation sequences
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
  else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
  else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
  else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}

  ## Reduced model: intercept plus all nuisance predictors
  Z <- X[, -tested, drop = FALSE]
  
  fit.reduced <- lme_fast(Y = Y,X = Z,id = random,add_intercept = FALSE)
  
  # Marginal residuals retain participant-level random effects.
  # Therefore permutations must preserve the participant covariance blocks.
  fitted0 <- Z %*% fit.reduced$coefficients
  residuals0 <- Y - fitted0
  
  shared <- mori::share(list(
    X = X,
    id = as.integer(random),
    fitted0 = fitted0,
    residuals0 = residuals0,
    permseq = permseq
  ))
  
  # shared remains referenced until all workers finish.
  rm(Z, fit.reduced, fitted0, residuals0,permseq, X, Y)

  cl <- parallel::makeCluster(
    min(nthread, nperm),
    type = "PSOCK"
  )
  
  pb <- NULL
  
  on.exit({
    if (!is.null(pb)) {
      try(close(pb), silent = TRUE)
    }
    
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  # Load mori before shared objects are deserialized in workers.
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
  
  start=Sys.time()
  message("\nEstimating permuted network strengths...\n")

  max.netstr <- foreach::foreach(
    perm = seq_len(nperm),
    .combine = "rbind",
    .packages = c("FCtools", "mori"),
    .options.snow = opts
  ) %dopar% {
    
    
    # Block-restricted residual permutation, shared across all edges.
    Yperm <- shared$fitted0 +shared$residuals0[shared$permseq[, perm],,drop = FALSE]
    fit.perm <- lme_fast(Y = Yperm,X = shared$X,id = shared$id,add_intercept = FALSE)
    
    # Do not silently discard failed permutations or treat them as zeros.
    if (anyNA(fit.perm$status) ||
        any(fit.perm$status != 0L)) {
      stop(
        "Nonzero fitting status in permutation ", perm,
        ". Inspect convergence before continuing."
      )
    }
    
    z.perm <- numeric(nedge)
    z.perm[active] <- as.numeric(
      fit.perm$t_stat[tested, ]
    )
    
    netstr <- cluster.stat(data = z.perm,nnodes = nnodes,tcrit = tcrit)
    
    if (is.null(dim(netstr))) {
      c(0, 0)
    } else {
      c(
        max(netstr[, 1L]),
        max(netstr[, 2L])
      )
    }
  }
  
  # foreach returns a vector when there is only one iteration.
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
      orig.clust[row, 3] <-round((1 + sum(max.netstr[, 1] >= orig.clust[row, 1])) /(nperm + 1),4)
      orig.clust[row, 4] <-round((1 + sum(max.netstr[, 2] >= orig.clust[row, 2])) /(nperm + 1),4)
    }

    #formatting results table
    orig.clust=cbind(c(1:nrow(orig.clust)),orig.clust)
    colnames(orig.clust)=c("network","strength.unweighted","strength.weighted","p.unweighted","p.weighted")
  
    #objects to return
    returnobj=list(orig.clust,z.orig, tcrit,max.netstr)
    names(returnobj)=c("results","t.orig","tcrit","max.netstr")
    return(returnobj)
  }
}


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
