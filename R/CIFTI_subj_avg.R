########################################################################################################
########################################################################################################
#' @title CIFTI_subj_avg
#'
#' @description Averaging multiple runs of vertex-wise(CIFTI) data within subject 
#'
#' @details Some fMRI datasets contain multiple runs per subject. Currently XCP-d generates ReHo or ALFF data for these multiple runs without concatenating the runs within a subject. This function averages the vertex-wise data such that each subject will only have one row of cifti data.
#'
#' @param rds The filepath of the `FSLRvextract`-generated .rds file
#' @param filename The filename of the output subject-averaged .rds file

#' @returns outputs a .rds file containing a matrix where each row represent each subject and each column represent each vertex.
#'
#' @examples
#' \dontrun{
#' CIFTI_subj_avg(rds = "COBRE_reho.rds", filename = "averaged/COBRE_reho.rds")
#' }
#' @export

########################################################################################################
########################################################################################################
CIFTI_subj_avg=function(rds, filename)
{
  dat=readRDS(rds)
  sublist.repeated=sub("\\_.*", "", dat[[1]])
  sublist=unique(sublist.repeated)
  sublist=sublist[order(sublist)]
  
  output.mat=matrix(NA, nrow=length(sublist), ncol=ncol(dat[[2]]))
  
  for(sub in 1:length(sublist))
  {
    idx=which(sublist.repeated==sublist[sub])
    if(length(idx)>1)
    {
      output.mat[sub,]=colMeans(dat[[2]][idx,])
    } else if(length(idx==1))
    {
      output.mat[sub,]=dat[[2]][idx,]
    }
  }
  print(sublist)
  saveRDS(output.mat,file = filename)
}
