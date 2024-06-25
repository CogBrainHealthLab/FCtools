#' @title extractFS.ROI
#'
#' @description function to extract Destrieux atlas and ASEG ROIs
#'
#' @details This function reads the *.table freesurfer output files and outputs two .csv files for the Destrieux and ASEG ROIs
#' @param basename The base name of the .csv file. Set to `dat` by default.
#' @returns Returns outputs two .csv files
#'}
#' @examples
#' \dontrun{
#' extractFS.ROI(basename="dat")
#' }
#' @export
##################################################################################################################
##################################################################################################################
##to train CPM models using a fixed p value threshold
extractFS.ROI=function(basename="dat")
{
  . <- LHDestrieux.table <- NULL 
  . <- RHDestrieux.table <- NULL 
  . <- aseg.vol.table <- NULL 
  lh=read.table("LHDestrieux.table",header = T)
  rh=read.table("RHDestrieux.table",header = T)
  CT=cbind(lh[,-c(1,76:78)],rh[,-c(1,76:78)])
  lh$eTIV
  subGM=read.table("aseg.vol.table",header = T)[,-c(1,15,34:40)]/lh$eTIV
  write.table(CT,file=paste0(basename,"_CT.csv"),sep =",",row.names = F)
  write.table(subGM,file=paste0(basename,"_subGM.csv"),sep =",",row.names = F)
}
