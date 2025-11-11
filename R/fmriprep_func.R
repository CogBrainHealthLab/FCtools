############################################################################################################################
############################################################################################################################
#' @title extract_headmotion
#'
#' @description Extracting head motion measurements.
#'
#' @details This function extracts head motion measurements from an fMRIprep output directory, and outputs these measurements in a .csv file
#'
#' @param filename the desired filename, with a *.csv extension, of the output.Set to 'motiondat.csv' by default
#' @param start calculate motion parameters from this frame index. set to 1 by default
#'
#' @returns outputs a .csv file containing columns of
#' \itemize{
#'  \item `subj` Subject ID.
#'  \item `no.frames` Number of frames.
#'  \item `FD` Framewise displacement.
#'  \item `RMSD` Root Mean Squared Displacement.
#'  \item `FD.20` Number of frames with FD>0.2.
#'  \item `RMSD.25` Number of frames with RMSD>0.25.
#'}
#' @examples
#' \dontrun{
#' extract_headmotion()
#' }
#' @export
########################################################################################################
########################################################################################################
extract_headmotion=function(filename="motiondat.csv",start=1)
{
  filelist = list.files(pattern = "_desc-confounds_timeseries.tsv",recursive = T)
  sub_task_list=gsub(pattern="_desc-confounds_timeseries.tsv",replacement = "",filelist)
  motiondat = data.frame(cbind(basename(sub_task_list), rep(0, NROW(filelist)),rep(0, NROW(filelist)), rep(0, NROW(filelist)), rep(0,NROW(filelist)), rep(0, NROW(filelist))))
  colnames(motiondat) = c("subj_task", "no.frames", "FD", "RMSD","FD.20", "RMSD.25")
  
  for (fileno in 1:NROW(filelist)) {
    cat(paste0(filelist[fileno],"\n"))
    confounds = read.table(file = filelist[fileno],sep = "\t", header = TRUE)
    if(start !=1){confounds=confounds[-c(1:(start-1)),]}
    confounds[confounds == "n/a"] = NA
    motiondat$RMSD[fileno] = mean(as.numeric(confounds$rmsd),na.rm = T)
    motiondat$no.frames = NROW(confounds)
    motiondat$FD[fileno] = mean(as.numeric(confounds$framewise_displacement),na.rm = T)
    motiondat$RMSD.25[fileno] = NROW(which(confounds$rmsd >0.25))/nrow(confounds)
    motiondat$FD.20[fileno] = NROW(which(confounds$framewise_displacement >0.2))/nrow(confounds)
  }
  write.table(motiondat, file = filename,sep = ",",quote = F,row.names = F)
  cat(paste0("headmotion data saved to ",filename))
}
