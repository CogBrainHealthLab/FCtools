#' @title aslprep.QC
#'
#' @description collating all aslprep QC data into a single .csv file
#'
#' @details This function collates all aslprep QC data into a single .csv file
#'
#' @param filename Filename of the .csv file. Set to `ASLPREP_QC.csv` by default
#' @returns outputs a .csv file containing the aslprep QC data
#'
#' @examples
#' \dontrun{
#' aslprep.QC()
#' }
#' @export

########################################################################################################
########################################################################################################

asl.QC=function(filename="ASLPREP_QC.csv")
{
  filelist=list.files(path = Sys.glob("sub-*"),pattern = "_desc-qualitycontrol_cbf.tsv", recursive = T,full.names = T)
  subses=basename(gsub(pattern = "_desc-qualitycontrol_cbf.tsv","",filelist))  
    for(sub in 1:length(filelist))
    {
      if(sub==1)
      {
      dat=read.table(filelist[sub],header = T)  
      } else
      {
        dat=rbind(dat,read.table(filelist[sub],header = T) )
      }
    }
  dat$subses=subses
  write.table(dat, file=filename, row.names = F, sep=",")
}
