############################################################################################################################
############################################################################################################################
#' @title extract_links
#'
#' @description Extracting links from data manifests
#'
#' @details This function extracts links from NDA formatted data manifest files
#'
#' @param manifest the filepath of the manifest file. Set to `datastructure_manifest.txt` by default.
#' @param files a vector of keywords that contained within the file names of the files to be downloaded.
#' @param filename The desired filename, with a *.txt extension, of the text file containing the links. Set to `download_list.txt` by default
#' @returns outputs a .txt file containing the links filtered from the data manifest file
#'
#' @examples
#' \dontrun{
#' extract_links(manifest="datastructure_manifest.txt", filename="download_list.txt", files=c("REST_Atlas_MSMAll_hp0_clean.dtseries.nii", "Movement_RelativeRMS.txt"))
#' }
#' @importFrom stringr str_detect
#' @export
########################################################################################################
########################################################################################################

extract_links=function(manifest="datastructure_manifest.txt",files,filename="downloadlist.txt")
{
  #read manifest and remove first row; the first row contains description of the column, hence not used.
  manifest=read.table(manifest,header = T)[-1,]
  
  #identify indices of the associated file
  idx.list=list()
  for(file in 1:NROW(files))
  {
    idx.list[[file]]=which(stringr::str_detect(pattern = files[file],string = manifest$associated_file)==T)
  }
  #check if files were found
  if(length(unique(unlist(idx.list)))==0)
  {
    stop("No files were found. Check your files argument")
  }
  #output filelist as a text file
  write.table(manifest$associated_file[unique(unlist(idx.list))],file=filename, quote = F, row.names = F, col.names = F)
}
########################################################################################################
########################################################################################################

#' @title extractFC
#'
#' @description Extracting FC matrices from HCP-preprocessed fMRI volumes
#'
#' @details This function extracts FC matrices from fMRI volumes preprocessed with the HCP pipeline with the option
#' of carrying out global signal regression and scrubbing. The directory structure and file names must adhere to the NDA/HCP data structure
#' and file naming conventions. \href{https://www.humanconnectome.org/software/get-connectome-workbench)}{Connectome Workbench} has to be
#' installed prior to running the function.
#'
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped
#' @param base_dir The filepath to directory containing the subject folders. Set to `fmriresults01/"` by default.
#' @param atlas The version (number of parcels) of the Schaefer atlas; it has to be in multiples of 100, from 100 to 1000. The specified atlas template will be automatically downloaded from here if it does not exist in the current directory.
#' @param FD option to use FD instead of RMSD for measuring head motion. Set to `FALSE` by default.
#' @param task The name of the task (case-sensitive), without any numbers or acquisition direction labels
#' @param extension The filename extension of the fMRI volumes. Set to `_Atlas_MSMAll_hp0_clean.dtseries.nii` by default
#' @param movement.extension The filename extension of the head motion files. Set to `Movement_RelativeRMS.txt` by default
#' @param motion.thresh the threshold (Root Mean Squared displacement) used for scrubbing frames. Set to 0.25 by default.
#' @param GSR If set to TRUE, global signal and its square-term will be regressed out from the fMRI timeseries data. Note, that while the use of GSR enhances brain-behavior association, its use is controversial in the literature
#' @param scrub If set to TRUE, remove frames with excessive head motion (relative RMS displacement> 0.25) and also remove segments in between two time-points if the segment has less than 5 frames. This is the ‘full scrubbing; method described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811913009117)}{Power et. al (2014)}. This is not recommended for task-based fMRI volumes.
#' @param output_dir The output directory where the FC matrices will be saved
#' @param report_filename The desired filename of the report containing columns of `subj`(subject ID),`Mean_RMS`(Mean Root Mean Squared displacement), and `no_frames_removed`(number of frames removed; only if `scrub=TRUE`)
#' @returns outputs M x M matrices in the `output_dir` and a report file (.csv format) in the working directory containing the head motion measurements
#'
#' @examples
#' \dontrun{
#' extractFC(wb_path = "/home/junhong.yu/workbench/",task="CONFLICT",atlas=100,GSR=TRUE,scrub=TRUE,output_dir = "FCmat")
#' }
#' @importFrom stringr str_detect
#' @importFrom ciftiTools ciftiTools.setOption read_cifti read_xifti newdata_xifti move_from_mwall
#' @importFrom fMRItools nuisance_regression
#' @export

########################################################################################################
########################################################################################################
extractFC=function(wb_path,
                   base_dir="fmriresults01/",
                   atlas,
                   FD=F,
                   motion.thresh=0.25,
                   scrub=F,
                   GSR=F,
                   task,
                   extension="_Atlas_MSMAll_hp0_clean.dtseries.nii",
                   movement.extension="Movement_RelativeRMS.txt",
                   output_dir="FCmat",
                   report_filename="report.csv")
{
  ##load and configure ciftitools
  
  ciftiTools::ciftiTools.setOption('wb_path', wb_path)
  
  ##create output_dir if it doesnt exist
  dir.create(file.path(output_dir), showWarnings = FALSE)
  
  ## file and subject listing
  fmri.filelist=list.files(path = base_dir,pattern = paste(task,".*",extension,sep=""),recursive = T,full.names=T)
  fmri.filelist=fmri.filelist[order(fmri.filelist)]
  sub.list=list.dirs(base_dir,recursive=F,full.names=F)
  movement.filelist=list.files(path = base_dir,pattern =movement.extension,recursive = T,full.names=T)
  movement.filelist=movement.filelist[order(movement.filelist)]
  
  ##setup report dataframe
  report=data.frame(matrix(NA,nrow=length(sub.list),ncol=2))
  if(FD==F)
  {
    colnames(report)=c("subj","Mean_RMS")
  } else
  {
    colnames(report)=c("subj","Mean_FD")
  }
  
  report$subj=sub.list
  if(scrub==T)
  {
    report$no_frames_removed=NA
  }
  ##checks
  #check missing fMRI volumes
  if(length(sub.list)>length(fmri.filelist))
  {
    missing.idx=rep(NA,length(sub.list))
    for(sub in 1:NROW(sub.list))
    {
      missing.idx[sub]=which(stringr::str_detect(pattern = sub.list[sub],string = fmri.filelist)==T)
    }
    cat("the number of fMRI volumes is less the number of subject folders. the following subjects do not have an fMRI volume\n")
    cat(paste(sub.list[which(is.na(c(missing.idx)))]))
    cat("\nYou can either delete these subject folders or re-download their files\n")
    stop()
  }
  
  #atlas check
  if(is.na(match(atlas,(1:10)*100)))
  {
    stop("Atlas should be a multiple of 100 from 100 to 1000")
  }
  #download atlas template if it is missing
  if(!file.exists(paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep="")))
  {
    cat(paste("The",paste("Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""),"template file does not exist and will be downloaded\n"),sep=" ")
    download.file(url=paste0("https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii"),
                  destfile =paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""),mode = "wb")
    
  }
  parc=c(as.matrix(read_cifti(paste(system.file(package='FCtools'),"/data/Schaefer2018_",atlas,"Parcels_7Networks_order.dlabel.nii",sep=""))))
  
  
  ## loop thru subjects
  for (sub in 1:NROW(sub.list))
  {
    start=Sys.time()
    cat(paste(sub.list[sub]))
    #filepaths
    fmri.path=fmri.filelist[which(stringr::str_detect(pattern = sub.list[sub],string = fmri.filelist)==T)]
    movement.path=movement.filelist[which(stringr::str_detect(pattern = sub.list[sub],string = movement.filelist)==T)]
    
    ##check number of movement files
    if(length(movement.path)>1)
    {
      movement.dat=list()
      frame.ends=list()
      for(volume in 1:length(movement.path))
      {
        if(FD==F)
        {
          movement.dat[[volume]]=read.table(movement.path[volume])$V1
        } else
        {
          movement.dat[[volume]]=extractFD(read.table(movement.path[volume])[,1:6])
        }
        
        if(volume!=length(movement.path))
        {
          frame.ends[[volume]]=c(NROW(movement.dat[[volume]]),NROW(movement.dat[[volume]])+1)
        }
      }
      for(volume in 2:(length(movement.path)-1))
      {
        frame.ends[[volume]]=frame.ends[[volume-1]][1]+frame.ends[[volume]]
      }
      movement.dat=unlist(movement.dat)
      frame.ends=unlist(frame.ends)
    } else
    {
      if(FD==F)
      {
        movement.dat=read.table(movement.path)$V1
      } else
      {
        movement.dat=extractFD(read.table(movement.path)[,1:6])
      }
    }
    if(FD==F)
    {
      report$Mean_RMS[sub]=mean(movement.dat)
    } else
    {
      report$Mean_FD[sub]=mean(movement.dat)
    }
    
    
    ##check number of fmri volumes, if multiple fmri volumes are detected, they are to be concatenated
    if(length(fmri.path)>1)
    {
      xii.list=list()
      for(vol in 1:length(fmri.path))
      {
        xii.list[[vol]]=ciftiTools::read_xifti(fmri.path[vol], brainstructures="all")
      }
      xii.base=xii.list[[1]]
      
      for(vol in 1:length(fmri.path))
      {
        if(vol==1)
        {
          xii.mat=as.matrix(xii.list[[1]])
        } else
        {
          xii.mat=cbind(xii.mat,as.matrix(xii.list[[vol]]))
        }
      }
      remove(xii.list,vol)
    } else
    {
      xii.base=ciftiTools::read_xifti(fmri.path, brainstructures="all")
      xii.mat=as.matrix(xii.base)
    }
    
    ##scrubbing
    if(scrub==T)
    {
      #identify frames for scrubbing
      exc_frames=which(movement.dat>motion.thresh)
      exc_frames=exc_frames[order(exc_frames)]
      exc_framesOLD=exc_frames
      
      for (frameno in 1:(NROW(exc_framesOLD)-1))
      {
        if(exc_framesOLD[frameno+1]-exc_framesOLD[frameno] <5 )
        {
          if(length(movement.path)>1)
          {
            totest=match(frame.ends,(exc_framesOLD[frameno]:exc_framesOLD[frameno+1]))
            totest=totest[-is.na(totest)]
            if(length(totest)!=2)
            {
              exc_frames=c(exc_frames,exc_framesOLD[frameno]:exc_framesOLD[frameno+1])
            }
          } else
          {
            exc_frames=c(exc_frames,exc_framesOLD[frameno]:exc_framesOLD[frameno+1])
          }
        }
      }
      exc_frames=unique(exc_frames)
      exc_frames=exc_frames[order(exc_frames)]
      
      #remove frames if necessary
      if(NCOL(xii.mat)!=length(movement.dat))
      {
        cat(paste("Number of timepoints in movement files do not match number of timepoints in fMRI volume.\nIt is likely that the subject is missing some movement files.\nScrubbing will not be carried out for",sub.list[sub]))
        xii.scrubbed=xii.mat
      } else
      {
        xii.scrubbed=xii.mat[,-exc_frames]
        report$no_frames_removed[sub]=length(exc_frames)
      }
      n_timepoints=NCOL(xii.scrubbed)
    } else
    {
      xii.scrubbed=xii.mat
      n_timepoints=NCOL(xii.scrubbed)
    }
    
    ##Global Signal Regression
    if(GSR==T)
    {
      dmat=xii.scrubbed
      gsig=colMeans(dmat)
      dmat=fMRItools::nuisance_regression(dmat, cbind(1, gsig,gsig^2))
      xii.regressed=ciftiTools::newdata_xifti(xii.base, dmat)
      xii.final=ciftiTools::move_from_mwall(xii.regressed, NA)
      remove(dmat,xii.base,xii.scrubbed,xii.regressed,gsig)
    } else if(GSR==F)
    {
      xii.final=ciftiTools::newdata_xifti(xii.base, xii.scrubbed)
      xii.final=ciftiTools::move_from_mwall(xii.final, NA)
      remove(xii.base,xii.scrubbed)
    }
    
    ##generate parcellated timeseries
    sub_keys=as.numeric(xii.final$meta$subcort$labels) - 2 +atlas
    brain_vec=c(parc, sub_keys)
    xii_pmean=matrix(ncol=atlas+19,nrow=n_timepoints)
    timeseries.dat=as.matrix(xii.final)
    
    for (p in 1:(atlas+19))
    {
      data_p=timeseries.dat[brain_vec==p,]
      xii_pmean[,p]=colMeans(data_p,na.rm = T)
    }
    #output FC matrix
    write.table(cor(xii_pmean), file=paste(output_dir,"/",sub.list[sub],".csv",sep=""), row.names = F, col.names = F, sep=",")
    remove(xii.final, sub_keys,brain_vec,xii_pmean,data_p)
    end=Sys.time()
    
    cat(paste(" completed in",round(difftime(end,start, units="secs"),1),"secs\n",sep=" "))
  }
  write.table(report,file = report_filename,sep = ",",row.names = F)
}
########################################################################################################
########################################################################################################
##function to extract FD, adapted from fMRItools

extractFD=function(mot_dat)
{
  detrend = F
  mot_dat[, 4:6] <- mot_dat[, 4:6]
  if (!isFALSE(detrend)) {
    if (isTRUE(detrend)) {
      detrend <- 4
    }
    mot_dat <- fMRItools::nuisance_regression(mot_dat, cbind(1, dct_bases(nrow(mot_dat), detrend)))
  }
  mot_datdiff <- apply(mot_dat, 2, diff, lag = 1)
  return(c(rep(0, 1), rowSums(abs(mot_datdiff))))
}
