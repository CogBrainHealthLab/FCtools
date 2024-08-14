## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title intersubject_similarity
#'
#' @description This function runs an intersubject similarity analysis to determine if there is a relationship between similarity in FC and similarity in one or multiple outcomes.
#'
#' @details This function runs an intersubject similarity analysis to determine if there is a relationship between similarity in FC and similarity in one or multiple outcomes. The outcome(s) will be z-standardized prior to calculating the intersubject similarity in the outcome(s).
#'
#' @param FC_data An N x E matrix containing the vectorized edges; where N = number of subjects, E=number of edges
#' @param outcome A numerical vector (single outcome) or matrix (multiple outcomes) containing the values of the outcome(s) of interest
#' @param mode When set to `"diff"`, FC similarity is calculated as absolute difference between the FC vectors of a pair of subjects. When set to `"corr"`,FC similarity is calculated as 1 - (pearson's correlation coefficient between the FC vectors of a pair of subjects). Set to `"diff"` by default.
#' @param nperm number of permutations for the correlation test between FC similarity and outcome similarity
#'
#' @returns a list object containing
#' \itemize{
#'  \item `FC_difference_matrix` The FC difference matrix
#'  \item `Outcome_difference_matrix` The outcome difference matrix
#'  \item `permutation_data` The permuted correlation values
#'  }
#' @examples
#'
#' results=runif(7021, min = -1, max = 1)
#' vizChord(data=results, filename="FC_chord119.png")
#'
#' \dontrun{
#' results=intersub(FC_data = dat_FC, outcome=dat_beh[,10:15],mode="diff")
#' }
#' @export
########################################################################################################
########################################################################################################
intersubject_similarity=function(FC_data, outcome,mode="diff", nperm=1000)
{
  outcome=data.matrix(outcome)
  
  #incomplete data check
  idxF=which(complete.cases(outcome)==F)
  if(length(idxF)>0)
  {
    cat(paste("outcome contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    outcome=outcome[-idxF,]
    FC_data=FC_data[-idxF,]
  }

  ##computing similarity matrices
  a=1
  n=NROW(FC_data)
  subjlist=1:n
  subjlist2=subjlist
  sim_dat=matrix(NA, ncol=4,nrow=((n*n)-n)/2)
  
  outcome=scale(outcome)
  for (subj1 in 1:(n-1))
  {
    subjlist2=subjlist2[-which(subjlist2==subj1)]
    
    if(mode=="diff")
    {
      for (subj2 in 1:(NROW(subjlist2)))
      {
        sim_dat[a,1]=subj1
        sim_dat[a,2]=subjlist2[subj2]
        sim_dat[a,3]=sum(abs(FC_data[subj1,]-FC_data[subjlist2[subj2],])) #similarity between FC
        sim_dat[a,4]=sum(abs(outcome[subj1,]-outcome[subjlist2[subj2],])) #similarity between behavioral outcomes
        a=a+1
      }  
    }
    else if(mode=="corr")
    {
      for (subj2 in 1:(NROW(subjlist2)))
      {
        sim_dat[a,1]=subj1
        sim_dat[a,2]=subjlist2[subj2]
        sim_dat[a,3]=1-cor(FC_data[subj1,],FC_data[subjlist2[subj2],]) #similarity between FC
        sim_dat[a,4]=sum(abs(outcome[subj1,]-outcome[subjlist2[subj2],])) #similarity between behavioral outcomes
        a=a+1
      } 
    }
  }
  ##computing permutations
  perm_dat=matrix(NA, ncol=nperm,nrow=((n*n)-n)/2)
   
  for(perm in 1:nperm)
  {
    a=1
    subjlist2=subjlist
    outcome.perm=data.matrix(outcome[sample(1:n),])
    for (subj1 in 1:(n-1))
    {
      subjlist2=subjlist2[-which(subjlist2==subj1)]
      
      for (subj2 in 1:(NROW(subjlist2)))
      {
        perm_dat[a,perm]=sum(abs(outcome.perm[subj1,]-outcome.perm[subjlist2[subj2],])) #similarity between behavioral outcomes
        a=a+1
      }  
    }  
  }
  perm.rho=as.numeric(cor(sim_dat[,4],perm_dat))
  rho=cor(sim_dat[,3],sim_dat[,4])
  
  ##format p values
  p=length(which(perm.rho>rho))/nperm
   if(p==0)
   {
     p=paste0("<",1/nperm)
   } else
   {
     p=paste0("=",p)
   }
  
  ##preparing objects to return
  FC.sim.mat=matrix(0,nrow=n,ncol=n)
  outcome.sim.mat=matrix(0,nrow=n,ncol=n)
  
  #converting upper triangle matrix to symmetric matrix
  FC.sim.mat[upper.tri(FC.sim.mat)]=sim_dat[,3]
  FC.sim.mat=FC.sim.mat+t(FC.sim.mat)
  FC.sim.mat[FC.sim.mat==0]=NA
  
  outcome.sim.mat[upper.tri(outcome.sim.mat)]=sim_dat[,4]
  outcome.sim.mat=outcome.sim.mat+t(outcome.sim.mat)
  outcome.sim.mat[outcome.sim.mat==0]=NA
  
  #output results as a message
  cat(paste0("\nCorrelation between FC and Outcome similarity matrice = ",round(rho,3)," ; p ",p))
  
  return.obj=list(FC.sim.mat,outcome.sim.mat,perm_dat)
  names(return.obj)=c("FC_difference_matrix","Outcome_difference_matrix","permutation_data")
  return(return.obj)
}
