
########################################################################################################
########################################################################################################

#' @title CIFTItoFC
#'
#' @description Extracting FC matrices from CIFTI files
#'
#' @details This function extracts FC matrices from CIFTI volumes postprocessed with XCP_d
#'
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped. Set to `/home/junhong.yu/workbench/bin_rh_linux64` by default
#' @param path The filepath to directory containing the subject folders. Set to `./.` by default.
#' @param atlas The version (number of parcels) of the Schaefer atlas; it has to be in multiples of 100, from 100 to 1000. set to `200` by default. The specified atlas template will be automatically downloaded from here if it does not exist in the current directory.
#' @param dtseries The filename extension of the fMRI volumes. Set to `_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii` by default
#' @param timeseries If set to `TRUE`, FC matrices will not be computed, instead the parcellated timeseries will be return in a list format, where each subject's parcellated timeseries data is contained within a list element. Set to `FALSE` by default.
#' @param concat_subj When set to `TRUE` (default), timeseries data from multiple runs/sessions of the same subject is first Z-scaled and then concatenated into a single larger timeseries data frame before computing the FC.
#' @param filename Filename of the concatenated FC vector file. Set to `FC.rds` by default
#' @returns outputs N (number of subjects) x E (number of unique) matrices as a .rds file 
#'
#' @examples
#' \dontrun{
#' CIFTItoFC(wb_path="bin_windows64/",dtseries="_ses-0._task-rest_acq-.._run-0._space-fsLR_den-91k_desc-denoised_bold.dtseries.nii",atlas = 200,filename="test.RDS")
#' }
#' @importFrom stringr str_detect

#@importFrom ciftiTools ciftiTools.setOption read_cifti read_xifti newdata_xifti move_from_mwall psych fisherz
#@export

########################################################################################################
########################################################################################################
CIFTItoFC=function(path="./",wb_path="/home/junhong.yu/workbench/bin_rh_linux64", dtseries="_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii", timeseries=F, concat_subj=TRUE, atlas=200,filename="FC.rds")
{
  filelist=list.files(path=path,pattern=dtseries,recursive = T)
  sublist=unique(gsub(pattern ="/ses-.*/func",replacement = "",dirname(filelist)))
  
  ##load and configure ciftitools
  ciftiTools::ciftiTools.setOption('wb_path', wb_path)
  
  ##atlas checks
  if(is.na(match(atlas,(1:10)*100)))  {stop("\nAtlas should be a multiple of 100 from 100 to 1000")}
  #download atlas template if it is missing
  if(!file.exists(paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep="")))
  {
    cat(paste("\nThe",paste("Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""),"template file does not exist and will be downloaded\n"),sep=" ")
    download.file(url=paste0("https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii"),
                  destfile =paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""),mode = "wb")
    
  }
  parc=c(as.matrix(ciftiTools::read_cifti(paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""),brainstructures=c("left","right"))))
  
  ## defining subcortical parcel indices for reordering ROIs subsequently
  reorder.subcortical.idx=c(9,18,8,17,6,3,13,1,11,10,19,7,16,5,15,4,14,2,12)+atlas
  
  if(concat_subj==F)  {sublist=filelist}
  
  if(timeseries==F)
  {
    all_FC=matrix(NA,nrow=length(sublist),ncol=(((atlas+19)^2)-(atlas+19))/2)  
  } else if(timeseries==T)
  {
    all_TS=list()
  }
  
  
  for (sub in 1:length(sublist))
  {
    start=Sys.time()
    cat(paste0("processing ",sublist[sub]," (",sub," / ",length(sublist),")..."))
    filelist.sub=filelist[stringr::str_detect(pattern=sublist[sub],string=filelist)]
    for (scan in 1:length(filelist.sub))
    {
      if(scan==1) 
      {
        xii0=ciftiTools::read_xifti(paste0(path,"/",filelist.sub[scan]), brainstructures="all")
        xii=scale(as.matrix(xii0))
        if(sum(colSums(xii)==0)!=0)
        {
          cat(paste0("WARNING: ",filelist.sub[scan]," contains column(s) with entirely 0 values\n"))
        }
      }
      else
      {
        xii=cbind(scale(as.matrix(ciftiTools::read_xifti(paste0(path,"/",filelist.sub[scan]), brainstructures="all"))),xii)
        if(sum(colSums(xii)==0)!=0)
        {
          cat(paste0("WARNING: ",filelist.sub[scan]," contains column(s) with entirely 0 values\n"))
        }
      }
    }
    
    if(length(filelist.sub)>1)
    {
      xii.final=ciftiTools::newdata_xifti(xii0, xii)
      timeseries.dat=as.matrix(ciftiTools::move_from_mwall(xii.final, NA))
    } else
    {
      xii.final=xii0
      timeseries.dat=as.matrix(ciftiTools::move_from_mwall(xii0, NA))
    }
    
    ##generate parcellated timeseries
    sub_keys=as.numeric(xii.final$meta$subcort$labels) - 2 +atlas
    brain_vec=c(parc, sub_keys)
    xii_pmean=matrix(ncol=atlas+19,nrow=ncol(xii))
    
    for (p in 1:(atlas+19))
    {
      xii_pmean[,p]=colMeans(timeseries.dat[which(brain_vec==p),],na.rm = T)
    }
    #reorder parcel indices for visualization purpose
    xii_pmean=xii_pmean[,c(1:atlas,reorder.subcortical.idx)] 
    
    if(timeseries==F)
    {
      FCmat=cor(xii_pmean)
      all_FC[sub,]=FCmat[upper.tri(FCmat,diag=F)]
      remove(FCmat)
    } else if(timeseries==T)
    {
      all_TS[[sub]]=xii_pmean
    }
  
    
    ##clean temp dir
    files_to_remove=list.files(tempdir(), full.names = TRUE)
    removedfiles=file.remove(files_to_remove)
    
    end=Sys.time()
    cat(paste0(" Completed in ", round(difftime(end,start, units="secs"),1),"s\n"))
    remove(xii, xii0, xii.final, xii_pmean, timeseries.dat,filelist.sub, start,end,removedfiles)
  }
  cat(paste0("Saving ", filename," ..."))
  
  if(timeseries==F)
  {
    if(concat_subj==T)
    {
      saveRDS(list(sublist,psych::fisherz(all_FC)),file=filename)  
    } else if (concat_subj==F)
    {
      saveRDS(list(basename(sublist),psych::fisherz(all_FC)),file=filename)  
    }
  } else if(timeseries==T)
  {
    if(concat_subj==T)
    {
      saveRDS(list(sublist,all_TS),file=filename)  
    } else if (concat_subj==F)
    {
      saveRDS(list(basename(sublist),all_TS),file=filename)  
    }
  }
}
