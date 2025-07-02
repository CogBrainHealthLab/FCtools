#' @title headmotion.XCP
#'
#' @description Calculating mean FD from XCP outputs and generating bash script to remove fMRI runs with excessive headmotion
#'
#' @details This function calculating mean FD from XCP outputs and generating bash script to remove fMRI runs with excessive headmotion. This is used mainly on the XCP-d post-processed HCP fMRI data.
#'
#' @param threshold mean FD threshold for identifying fMRI runs for deletion due to excessive headmotion. Set to `0.25` by default
#' @param filename Filename of the output.csv file. Set to `FD.csv` by default
#' @returns outputs a .csv file containing the columns of fMRI runs and FD, and a bash script `del.sh`
#'
#' @examples
#' \dontrun{
#' headmotion.XCP()
#' }
#' @export

########################################################################################################
########################################################################################################

headmotion.XCP=function(filename="FD.csv", threshold=0.25)
{
  filelist=list.files(pattern="_motion.tsv", recursive=T)
  FD.all=data.frame(matrix(NA,nrow=length(filelist),ncol=2))
  colnames(FD.all)=c("fMRI_run","FD")
  
  for (file in 1:length(filelist))
  {
    FD.all$fMRI_run[file]=basename(filelist[file])
    FD.all$FD[file]=mean(read.table(filelist[file],header=T, sep="\t")$framewise_displacement)
  }
  
  files.del=FD.all[which(FD.all$FD>threshold),1]
  del.sh=paste0("rm -rf ",gsub("_motion.tsv","*",files.del))

  if(length(files.del))>1)
  {
    cat(paste0(length(files.del), " fMRI runs had mean FD values >", threshold)) 
    write.table(del.sh,row.names=F, col.names=F,quote=F, file="del.sh")
  } else
  { cat(paste0("None of the fMRI runs had mean FD values >",threshold)) }
  FD.all$fMRI_run=gsub(pattern="_motion.tsv",replacement="",x=FD.all$fMRI_run)
  write.table(FD.all, row.names=F, quote=F, sep=",", file=filename)
}

