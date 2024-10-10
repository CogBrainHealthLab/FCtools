#' @title cpm.train
#'
#' @description Training a connectome-based predictive model using a fixed p-value threshold
#'
#' @details This function implements the connectome-based prediction model described in \href{https://www.nature.com/articles/nprot.2016.178)}{Shen et al. (2017)}
#' @param data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param outcome The outcome variable to predict
#' @param p The p-value threshold of a Pearson's correlation test between the feature and outcome that determines if the feature is selected. Set to 0.05 by default.
#' @returns Returns a list object containing
#' \itemize{
#'   \item `weights` vector of 1s (positive edges), 0s and -1s (negative edges) indicating the selected edges
#'   \item `pos.network.coef` linear regression coefficients of the positive network strength model
#'   \item `neg.network.coef` linear regression coefficients of the negative network strength model
#'   \item `both.network.coef` linear regression coefficients of the positive and negative network strength model
#'}
#' @examples
#' \dontrun{
#' model1=cpm.train(data=FC_data,outcome=dat_beh$age, p=0.05)
#' }
#' @export
##################################################################################################################
##################################################################################################################
##to train CPM models using a fixed p value threshold

cpm.train=function(data,outcome,p=0.05)
{
  ##checks
  #pvalues
  data=data.matrix(data)
  if(length(p)==1)  {p.posneg=c(p,p)}
  else if (length(p)==2)  {p.posneg=p}
  else {stop("Only one or two pvalues should be entered")}

  #number of rows match
  if(NROW(data)!=NROW(outcome))  {stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))}

  #missing data
  idx.missing=which(is.na(outcome)==T)
  if(NROW(idx.missing)>0)
  {
    cat(paste("\n",NROW(idx.missing)," missing values are detected in the outcome variable. Subjects with missing values will be excluded in the training procedure\n",sep=""))
    data=data[-idx.missing,]
    outcome=outcome[-idx.missing]
  }

  ##feature selection
  #critical values
  pos.tcrit=qt((p.posneg[1]/2), NROW(outcome)-2, lower=FALSE)
  neg.tcrit=qt((p.posneg[2]/2), NROW(outcome)-2, lower=FALSE)
  pos.rcrit=sqrt(pos.tcrit^2/(pos.tcrit^2+NROW(outcome)-2))
  neg.rcrit=sqrt(neg.tcrit^2/(neg.tcrit^2+NROW(outcome)-2))

  #binarizing pearsons r values
  r.mat=cor(data,outcome)
  weights=rep(0,NCOL(data))
  weights[r.mat> pos.rcrit]=1
  weights[r.mat< -neg.rcrit]=-1

  ##network models
  #positive model
  if(NROW(which(weights==1))>1) ## proceed only if at least 2 edges are selected
  {
    pos.netstr=rowSums(data[,which(weights==1)])
    pos.netstr.coef=lm(outcome~pos.netstr)$coefficients
  } else
  {
    pos.netstr.coef=NA
    cat("\nNone of edges are significantly and positively correlated with the outcome. The positive network model cannot be constructed\n")
  }

  #negative model
  if(NROW(which(weights==-1))>1) ## proceed only if at least 2 edges are selected
  {
    neg.netstr=rowSums(data[,which(weights==-1)])
    neg.netstr.coef=lm(outcome~neg.netstr)$coefficients
  } else
  {
    neg.netstr.coef=NA
    cat("\nNone of edges are significantly and negatively correlated with the outcome. The negative network model cannot be constructed\n")
  }

  # positive + negative model
  if(NROW(which(weights==-1))>1 & NROW(which(weights==1))>1) ## proceed only if at least 2 edges are selected in each of the earlier models
  {both.netstr.coef=lm(outcome~pos.netstr+neg.netstr)$coefficients}
  else {both.netstr.coef=NA}

  ##listing objects to return
  model=list(weights,pos.netstr.coef,neg.netstr.coef,both.netstr.coef)
  names(model)=c("weights","pos.network.coef","neg.network.coef","both.network.coef")
  return(model)
}
##################################################################################################################
##################################################################################################################
#' @title cpm.predict
#'
#' @description predicting scores using a previously trained CPM model
#'
#' @details This function takes a previously trained CPM model's weights to calculate the positive and negative network strengths for each subject.
#' Then these network strengths are multiplied with the regression coefficients in the CPM model to compute a predicted score.
#' @param model A list object generated using `cpm.train()`
#' @param test.data An N x E matrix containing the vectorized edges of subjects whose outcomes/scores are to be predicted. This matrix has to contain the same number of columns as the matrix of features used to train the CPM model
#' @param network A string to specify which network strength(`'positive'`,`'negative'` or `'both'`) to use in the outcome prediction. Set to `'both'` by default
#' @returns Returns a vector of predicted outcomes if `network='positive'` or `network='negative'`.
#' In the case of `network='both'`, a data.frame containing the predicted scores using the positive, negative, and the combined positive + negative models will be returned
#' @examples
#' \dontrun{
#' model1=cpm.predict(model=model1, data=FC_data.test,network="both")
#' }
#' @export

##################################################################################################################
##################################################################################################################
##to predict scores from previously generated cpm.train() models

cpm.predict=function(model,test.data, network="both")
{
  test.data=data.matrix(test.data)
  ##checks
  #number of rows match
  if(NROW(model$weights)!=NCOL(test.data))
  {stop(paste("\nThe number of predictors in the training data (",NROW(model$weights),") and testing data(",NCOL(test.data),")do not match!\n",sep=""))}

  ##select model {compute predscore}
  if(network=="positive")  {predscore=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]}
  else if(network=="negative")  {predscore=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]}
  else if(network=="both")
  {
    #check if positive and negative models are valid, and proceed accordingly
    if(is.na(model$pos.network.coef)[1])  {positive=rep(NA,NROW(test.data))}
    else {positive=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]}

    if(is.na(model$neg.network.coef)[1])  {negative=rep(NA,NROW(test.data))}
    else {negative=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]}

    if(is.na(model$pos.network.coef)[1] | is.na(model$neg.network.coef)[1])  {both=rep(NA,NROW(test.data))}
    else  {both=rowSums(test.data[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(test.data[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]}

    predscore=data.frame(positive,negative,both)
  }
  return(predscore)
}
##################################################################################################################
##################################################################################################################
#' @title cpm.train.cv
#'
#' @description K-fold cross-validation procedure to determine the optimal p-value threshold for feature selection
#'
#' @details This function runs the `cpm.train()` using a range of p-values and determines the optimal p-values that maximize predicted-actual correlation in a K-fold cross-validation paradigm.
#' @param data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param outcome The outcome variable to predict
#' @param p a vector of p-values to tested. If `p` is not specified, an automatically determined geometric sequence of p-values will be tested
#' @param nfolds number of cross-validation folds. Set to 5 by default
#' @returns Returns a list object containing
#' \itemize{
#'  \item `opt.pvals` A vector containing the 2 optimal p-values, one each for the positive and negative network models
#'  \item `results` A P(number of p-values tested) x 2(positive and negative network models) matrix containing the predicted-actual correlations for each of the p-value selection thresholds in each network model.
#'  \item `pvals` A vector of the p-value threshold tested
#'  }
#' @examples
#' \dontrun{
#' model1.cv=cpm.train.cv(data=FC_data,outcome=dat_beh$age)
#' }
#' @importFrom caret createFolds
#' @export

cpm.train.cv=function(data,outcome,p,nfolds=5)
{
  ##checks
  #number of rows match
  if(NROW(data)!=NROW(outcome))  {stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))}

  #missing data
  idx.missing=which(is.na(outcome)==T)
  if(NROW(idx.missing)>0)
  {
    data=data[-idx.missing,]
    outcome=outcome[-idx.missing]
  }
  ##generate geometric sequence of p-values from distribution of p values in the data
  if(missing(p))
  {
    #converting r to p
    r_to_p=function(r)
    {
      t=(r*sqrt(NROW(outcome)-2))/(sqrt(1-(r^2)))
      return(2 * (1 - pt(abs(t), NROW(outcome)-2)))
    }

    r.thresh=0.15 ##lower cutoff
    r.mat=cor(data,outcome)

    #separating p values for positive and negative rs
    pos.r.mat=r.mat[r.mat>r.thresh]
    neg.r.mat=r.mat[r.mat< (-r.thresh)]

    ##generating p values
    #order distribution of p values in increasing magnitude
    pos.r.mat=pos.r.mat[order(pos.r.mat)]
    neg.r.mat=neg.r.mat[order(-neg.r.mat)]

    #divide p values across 47 intervals
    n.int=47
    pos.interval=NROW(pos.r.mat)/n.int
    neg.interval=NROW(neg.r.mat)/n.int

    #select intervals from a geometric step progression; steps become progressively smaller
    intervals=c(9,17,24,30,35,39,42,44,45)
    p=matrix(NA,nrow=length(intervals)+1, ncol=2)

    #first pair of p values set according to r=+/-0.15 (lower boundary)
    p[1,]=c(r_to_p(r.thresh),r_to_p(-r.thresh))

    #iterating p values across geometric stepped intervals for positive and negative models
    p[2:NROW(p),1]=r_to_p(pos.r.mat[round(pos.interval*intervals)])
    p[2:NROW(p),2]=r_to_p(neg.r.mat[round(neg.interval*intervals)])

  } else
  {
    #checks for user-defined p values
    if(length(p)==1)  {stop("At least 2 p-values should be entered")}
    else if (NCOL(p)==1)  {p=cbind(p,p)}
  }

  ##setup CV folds
  folds=caret::createFolds(outcome,k=nfolds)

  ##training

  rvals=matrix(NA,nrow=NROW(p),ncol=2)
  for(iter in 1:NROW(p))
  {
    p.iter=p[iter,]
    for (fold in 1:length(folds))
    {
      #leaving a fold out
      train_outcome=outcome[-folds[[fold]]]
      train_dat=data[-folds[[fold]],]

      #training the CPM and applying model to predict scores
      train.mod=cpm.train(data = train_dat,outcome = train_outcome,p=p.iter)
      predscores.fold=cpm.predict(train.mod, test.dat=data[folds[[fold]],])

      #saving the scores from each fold to a single vector
      if(fold==1)  {predscores.all=predscores.fold}
      else  {predscores.all=rbind(predscores.all,predscores.fold)}
    }
    #check if positive and negative models are valid, and proceed accordingly
    if(anyNA(predscores.all[,1]))  {pos.net.r=NA}
    else {pos.net.r=cor(outcome[unlist(folds)],predscores.all[,1])}

    if(anyNA(predscores.all[,2]))  {neg.net.r=NA}
    else {neg.net.r=cor(outcome[unlist(folds)],predscores.all[,2])}

    #saving the r values for the iteration
    rvals[iter,]=c(pos.net.r,neg.net.r)
  }

  #inform user if p values are too small such that no edges are selected; only if user-defined p values are entered
  idx.NA.pos=which(is.na(rvals[,1]))
  if(length(idx.NA.pos)>0)
  {cat(paste("\nNote: No positive edges were selected when the p value of ",p[min(idx.NA.pos)]," or smaller was used", sep=""))}

  idx.NA.neg=which(is.na(rvals[,2]))
  if(length(idx.NA.neg)>0)
  {cat(paste("\nNote: No negative edges were selected when the p value of ",p[min(idx.NA.neg)]," or smaller was used", sep=""))}

  #identify the indices for p-values that produce the most accurate predictions
  r.pos.min.idx=which(rvals[,1]==max(rvals[,1],na.rm = T))
  r.neg.min.idx=which(rvals[,2]==max(rvals[,2],na.rm = T))

  #listing out objects to return
  results=list(c(p[r.pos.min.idx,1],p[r.neg.min.idx,2]),rvals,p)
  names(results)=c("opt.pvals","results","pvals")

  return(results)
}
##################################################################################################################
##################################################################################################################

#' @title cpm.lesion
#'
#' @description CPM with the leave-one-network-out lesion approach
#'
#' @details This function runs the `cpm.train()` and `cpm.predict()` while leaving out a network of edges each time.
#' The changes in predicted-actual correlation while a network being left out can be interpreted to be an
#' indication of the contribution of the network to predicting the outcome
#' @param train.data An N x E matrix containing the vectorized edges in the training dataset; where N = number of subjects, E=number of edges
#' @param test.data An N x E matrix containing the vectorized edges in the training dataset
#' @param train.outcome The outcome variable to predict, within the training dataset
#' @param test.outcome The outcome variable to predict, within the testing dataset
#' @param p The p-value threshold of a Pearson's correlation test between the feature and outcome that determines if the feature is selected. Set to 0.05 by default.
#' @returns Returns a matrix with the following columns
#' \itemize{
#'  \item `lesion.model` The CPM model in which the listed network is left out
#'  \item `positive` The predicted-actual correlations for each of the `lesion.model`s' positive network model
#'  \item `negative` The predicted-actual correlations for each of the `lesion.model`s' negative network model
#'  \item `both` The predicted-actual correlations for each of the `lesion.model`s' combined positive + negative network model
#'  }
#' @examples
#' \dontrun{
#' model1.lesion=cpm.lesion(train.data=FC_data.train,test.data=FC_data.test,train.outcome=train_dat$age, test.outcome=test_dat$age,p=0.05)
#' }
#'
#' @importFrom caret createFolds
#' @export
##Using the lesion approach (leave-one-network out) for CPM
cpm.lesion=function(train.data,test.data,train.outcome, test.outcome,p=0.05)
{
  data=data.matrix(dat_FC)

  ##atlas selection
  edge_lengths=c(4005,7021,23871,30135)
  if(is.na(match(NCOL(train.data),edge_lengths)))
  {
    stop("The length of the input vector does not fit any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")
  } else
  {
    atlas=match(NCOL(train.data),edge_lengths)
  }

  ##preparing atlas labels for removing networks of edges
  labels_dat=get('labels_dat')
  label=labels_dat[[match(NCOL(train.data),edge_lengths)]]
  nnode=NROW(label)
  networks.list=data.frame(unique(cbind(as.numeric(label$region),label$regionlabel)))
  names(networks.list)=c("netno","network.name")
  networks.list=networks.list[order(networks.list$netno),]

  ##CPM
  results=matrix(NA,nrow=NROW(networks.list)+1,ncol=4)

  #training the CPM (without any network exclusions) and applying model to predict scores
  model.allnetworks=cpm.train(data=train.data, outcome=train.outcome, p=p)
  pred.allnetwork=cpm.predict(model = model.allnetworks, test.data=test.data)
  results[1,2:4]=cor(test.outcome,pred.allnetwork)

  #CPM with one network removed each time
  for (netno in 1:NROW(networks.list))
  {
    #identifying indice of edges to remove
    FC_matrix=array(rep(NA,nnode^2),dim=c(nnode,nnode))
    FC_matrix[upper.tri(FC_matrix, diag=FALSE)] = 1:edge_lengths[atlas]
    FC_matrix.1net=FC_matrix
    FC_matrix.1net[which(label$region==netno),which(label$region==netno)]=NA
    edge.column=FC_matrix.1net[upper.tri(FC_matrix.1net, diag=FALSE)]
    remove.idx=which(is.na(edge.column)==T)

    #training the CPM and applying model to predict scores
    model.1net=cpm.train(data=train.data[,-remove.idx], outcome=train.outcome, p=p)
    pred.1net=cpm.predict(model = model.1net, test.data=test.data[,-remove.idx])
    results[netno+1,2:4]=cor(test.outcome,pred.1net)
  }

  ##saving results in a data.frame object for returning
  results[1,1]="none removed"
  results[c(2:(NROW(networks.list)+1)),1]=paste("removed",networks.list$network.name, sep=" ")
  results=data.frame(results)
  results[, c(2:4)]=sapply(results[, c(2:4)], as.numeric)
  names(results)=c("lesion.model","positive","negative","both")
  return(results)
}
##################################################################################################################
##################################################################################################################
