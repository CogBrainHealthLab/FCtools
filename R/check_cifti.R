########################################################################################################
########################################################################################################
#' @title check_cifti
#'
#' @description function to count number of 0 columns in cifti files
#'
#' @details This function searches for all *_space-fsLR_den-91k_bold.dtseries.nii files recursively. Then these files are read one at a time, data from the left and right cortices and the subcortex are then `rbind` to form a 91282 x N matrix, where N is the number of frames. This matrix is summed rowwise, and the number of rowwise sums that are equal to 0 are counted. If there are no issues with the structural-functional alignment, there should only be 546 zero values
#'
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped. Set to `/home/junhong.yu/workbench/bin_rh_linux64` by default
#' @param filename check the directory structure to ensure that each fMRI directory should only contain one fMRI file and at least one movement file. Set to `TRUE` by default.
#' @returns outputs a .csv file with a N x 2 data.frame, with the columns of `filename` and `Non-zero_values``
#'
#' @examples
#' \dontrun{
#' check_cifti()
#' }
#' @importFrom ciftiTools ciftiTools.setOption read_xifti
#' @export

########################################################################################################
########################################################################################################
check_cifti=function(filename="check_cifti.csv", wb_path="/home/junhong.yu/workbench/bin_rh_linux64")
{
  ciftiTools::ciftiTools.setOption('wb_path', wb_path)

  cat("Searching for *_space-fsLR_den-91k_bold.dtseries.nii files...\n")
  dirs=list.dirs(recursive = F)
  dirs=grep("sub-*", dirs, value = TRUE)
  dirs=gsub(pattern = "./", replacement="", x=dirs)

  for(sub in 1:length(dirs))
  {
    fmri.filelist.sub=list.files(path = dirs[sub],pattern = "*_space-fsLR_den-91k_bold.dtseries.nii",recursive = T,full.names=T)
    if(length(fmri.filelist.sub)>=1)
    {
      if(sub==1)    
      {
        fmri.filelist=fmri.filelist.sub
      } else
      {
        fmri.filelist=c(fmri.filelist,fmri.filelist.sub)
      }
    } else
    {
      warning(paste0(dirs[sub], "does not contain any *_space-fsLR_den-91k_bold.dtseries.nii files \n"))
    }
  }

  cat("Search completed\n")
  check=matrix(NA, ncol=2,nrow=length(fmri.filelist))
  for(file in 1:length(fmri.filelist))
  {
    xii=ciftiTools::read_xifti(fmri.filelist[1], brainstructures="all")  
    xii.combined=rbind(xii$data$cortex_left,xii$data$cortex_right,xii$data$subcort)
    check[file,1]=basename(fmri.filelist[file])
    rowsum=rowSums(xii.combined)
    check[file,2]=91282-length(which(abs(rowsum)>0))
    cat(paste0(fmri.filelist[file]," ",check[file,2], " Non-zero values\n"))
    remove(xii, xii.combined)
  }
  check=data.frame(check)
  colnames(check)=c("filename","Non-zero_values")
  write.table(check,file=filename,row.names = F, sep=",")
}

