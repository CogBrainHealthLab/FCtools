############################################################################################################################
############################################################################################################################
#' @title extract_linksXCP
#'
#' @description Extracting links of XCP-d minimal inputs and FreeSurfer surface data from data manifests
#'
#' @details This function extracts links of XCP-d minimal inputs files and FreeSurfer surface data from NDA formatted data manifest files
#'
#' @param manifest the filepath of the manifest file. Set to `datastructure_manifest.txt` by default.
#' @param task name (case-sensitive) of the task, without numbers or directional labels. E.g., `REST`
#' @param surf name (case-sensitive) of the FreeSurfer surface data. E.g., `thickness`
#' @param filename The desired filename, with a *.txt extension, of the text file containing the links. Set to `download_list.txt` by default
#' @param subjects If specified the only download links from these subjects will be included
#' @returns outputs a .txt file containing the links filtered from the data manifest file
#'
#' @examples
#' \dontrun{
#' extract_linksXCP(manifest = "datastructure_manifest.txt", task="REST",freesurfer = "thickness",subjects = c("0891566","0571649"))
#' }
#' @importFrom stringr str_detect
#' @export
########################################################################################################
########################################################################################################


extract_linksXCP=function(manifest="datastructure_manifest.txt",task,surf,filename="downloadlist.txt", subjects)
{
  files=list()
  if(!missing(task))
  {
    files=c("Atlas_MSMAll.dtseries.nii",
            "brainmask_fs.2.nii.gz",
            "_SBRef.nii.gz",
            "Movement_Regressors.txt",
            "Movement_AbsoluteRMS.txt",
            "MNINonLinear/T1w.nii.gz",
            "MNINonLinear/ribbon.nii.gz",
            "MNINonLinear/aparc\\+aseg",
            "MNINonLinear/T1w.nii.gz",
            "MNINonLinear/brainmask_fs.nii.gz",
            paste0("_",task,"?_...nii.gz"))
  }
  if(!missing(surf) & !missing(task))
  {
    files=c(files,"lh.cortex.label",
            "rh.cortex.label",
            "lh.sphere.reg",
            "rh.sphere.reg",
            paste0("lh.",surf),
            paste0("rh.",surf),
            "aseg.stats")
  } else if(!missing(surf) & missing(task))
  {
    files=c("lh.cortex.label",
            "rh.cortex.label",
            "lh.sphere.reg",
            "rh.sphere.reg",
            paste0("lh.",surf),
            paste0("rh.",surf),
            "aseg.stats")
  }
  #read manifest and remove first row; the first row contains description of the column, hence not used.
  manifest=read.table(manifest,header = T)[-1,]
  
  #identify indices of the associated file
  idx.list=list()
  for(file in 1:NROW(files))
  {
    idx.list[[file]]=which(stringr::str_detect(pattern = files[file],string = manifest$associated_file)==T)
  }
  
  filelist=manifest$associated_file[unique(unlist(idx.list))]
  
  #select subjects
  if(!missing(subjects))
  {
    for(sub in 1:length(subjects))
    {
      if(sub==1)
      {
        filelist.sel=filelist[stringr::str_detect(string=filelist,pattern =subjects[sub])]
      } else
      {
        filelist.sel=c(filelist.sel,filelist[stringr::str_detect(string=filelist,pattern =subjects[sub])])  
      }
    }
  } else
  {
    filelist.sel.alltasks=filelist
  }
  
  #remove MNINonLinear/brainmask_fs.2.nii.gz
  if(!missing(task))
  {
    filelist.sel.alltasks=filelist.sel.alltasks[-stringr::str_detect(string=filelist.sel.alltasks,pattern ="MNINonLinear/brainmask_fs.2.nii.gz")]    
    #filter tasks
    for (t in 1:length(task))
    {
    if(t==1)
      {
        filelist.sel=filelist.sel.alltasks[stringr::str_detect(string=filelist.sel.alltasks,pattern =task[1])]      
      } else
      {
        filelist.sel=c(filelist.sel,filelist.sel.alltasks[stringr::str_detect(string=filelist.sel.alltasks,pattern =task[t])])      
      }
    }
  }
  #check if files were found
  if(length(filelist.sel)==0)
  {
    stop("No files were found, check your task name")
  }
  #output filelist as a text file
  write.table(filelist.sel,file=filename, quote = F, row.names = F, col.names = F)
}


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

##currently this script concatenates the timeseries across multiple runs before denoising. 
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
#' @param task The name of the task (case-sensitive), without any numbers or acquisition direction labels
#' @param extension The filename extension of the fMRI volumes. Set to `_Atlas_MSMAll_hp0_clean.dtseries.nii` by default
#' @param movement.extension The filename extension of the head motion files. Set to `Movement_Regressors.txt` by default. Note the head motion file should be a text file containing either a single column of numbers (RMS measurements), or at least 6 columns (3 displacement + 3 rotational vectors). If there are more than 6 columns, only the first 6 columns will be used to compute framewise displacement.
#' @param motion.thresh the threshold used for scrubbing frames. Set to 0.25 by default.
#' @param motion.reg When set to `TRUE`, the 6 motion parameters (3 translation and 3 rotation), their derivatives and the squares of these parameters and their derivatives will be regressed out of the fMRI timeseries data. This motion regression procedure will only work if each subject folder contains a `Movement_Regressors.txt` file.
#' @param GSR If set to TRUE, global signal, its derivative and the squares of the global signal and its derivatives will be regressed out from the fMRI timeseries data.
#' @param bandpass Lower and upper limits of band pass filter. Signals below the lower limit and above the upper limit of the band pass filter will be eliminated. set to `c(0.01,0.1)` by default. Can turned off by setting to `NA`.
#' @param TR Temporal resolution of the fMRI data. Only used if bandpass is not `NA`.
#' @param scrub If set to TRUE, remove frames with excessive head motion (FD > 0.25) and also remove segments in between two time-points if the segment has less than 5 frames. This is the ‘full scrubbing; method described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811913009117)}{Power et. al (2014)}. This is not recommended for task-based fMRI volumes.
#' @param output_dir The output directory where the FC matrices and headmotion parameters will be saved. A default `output_dir` directory will be created if it does not exist.
#' @param output.timeseries When set to `TRUE`, the parcellated timeseries data will be output instead of the FC matrix.
#' @param overwrite If set to `FALSE`, subjects that already have their FC matrices in `output_dir` will be skipped. Set to `TRUE` by default.
#' @param dir.check check the directory structure to ensure that each fMRI directory should only contain one fMRI file and at least one movement file. Set to `TRUE` by default.
#' @returns outputs M x M matrices in `output_dir\FCmat`, headmotion measurements in `output_dir\headmotion` and a summary headmotion measurements across all subjects `output_dir\report.csv`.
#'
#' @examples
#' \dontrun{
#' extractFC(wb_path = "/home/junhong.yu/workbench/",task="REST",atlas=100,GSR=TRUE,scrub=TRUE,output_dir = "FCmat")
#' }
#' @importFrom stringr str_detect
#' @importFrom ciftiTools ciftiTools.setOption read_cifti read_xifti newdata_xifti move_from_mwall
#' @importFrom fMRItools nuisance_regression dct_bases
#' @export

########################################################################################################
########################################################################################################
extractFC=function(wb_path,
                   base_dir="fmriresults01/",
                   atlas,
                   subjects,
                   motion.thresh=0.25,
                   scrub=F,
                   motion.reg=T,
                   GSR=F,
                   bandpass=c(0.01,0.1),
                   TR,
                   task,
                   extension="_Atlas_MSMAll_hp0_clean.dtseries.nii",
                   movement.extension="Movement_Regressors.txt",
                   output_dir="output_dir",
                   output.timeseries=F,
                   overwrite=TRUE,
                   dir.check=TRUE)
{
  ##check base_dir and sub.list
  
  cat("ExtractFC tool: last updated 15 Jan 2025 6.00pm\n checking directory structure...\n\n")
  if(!dir.exists(base_dir))  {stop(paste("The base directory '",base_dir,"' does not exist. Please check if you are at the correct working directory",sep=""))}
  sub.list=list.dirs(base_dir,recursive=F,full.names=F)
  N.orig=length(sub.list)
  if(length(sub.list)==0) {stop("no subject folders were found. Please check if the 'base_dir' is correctly specified")}
  sub.list=sub.list[order(sub.list)]
  
  ##check bandpass filter
  if(!is.na(bandpass) & length(bandpass)!=2)  {stop("Band pass filter requires lower and upper limits of the frequency band to be specified")} 
  if(length(bandpass)==2 & missing(TR)) {stop("TR needs to be specified in order to apply band pass filtering")} 
  
  ##create output_dir and headmotion_dir if they dont exist
  dir.create(file.path(output_dir), showWarnings = FALSE)
  setwd(output_dir)
  dir.create(file.path("FCmat"), showWarnings = FALSE)
  dir.create(file.path("headmotion"), showWarnings = FALSE)
  setwd("..")
  
  ## file and subject listing
  fmri.filelist=list.files(path = base_dir,pattern = paste(task,".*",extension,sep=""),recursive = T,full.names=T)
  movement.filelist=list.files(path = base_dir,pattern =movement.extension,recursive = T,full.names=T)
  dir.list=list.dirs(base_dir,recursive=T,full.names=T)
  
  #if subjects vector is specified
  if(!missing(subjects))
  {
    fmri.path.idx=list()
    movement.path.idx=list()
    sub.list.idx=list()
    dir.list.idx=list()
    
    for (subj in 1:length(subjects))
    {
      fmri.path.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = fmri.filelist)==T)
      movement.path.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = movement.filelist)==T)
      sub.list.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = sub.list)==T)
      dir.list.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = dir.list)==T)
    }
    fmri.path.idx=unlist(fmri.path.idx)
    movement.path.idx=unlist(movement.path.idx)
    sub.list.idx=unlist(sub.list.idx)
    dir.list.idx=unlist(dir.list.idx)
    
    #if the user provides a wrong subject ID
    if(length(sub.list.idx)!=length(subjects))  {stop("not all specified subjects' directories could be found")}
    
    fmri.filelist=fmri.filelist[fmri.path.idx]
    movement.filelist=movement.filelist[movement.path.idx]
    dir.list=dir.list[dir.list.idx]
    sub.list=subjects
  }
  
  ##check if any fMRI volumes/movement files are found; returns error if 0 files are found
  if(length(fmri.filelist)==0)  {stop("no fMRI volumes were found. Please check if 'task' and 'extension' are correctly specified")}
  if(length(movement.filelist)==0) {stop("no movement files were found. Please check if 'movement.extension' is correctly specified")}
  
  fmri.filelist=fmri.filelist[order(fmri.filelist)]
  movement.filelist=movement.filelist[order(movement.filelist)]
  
  ##file size checks
  filesizes=file.info(fmri.filelist)$size
  if(any(filesizes<10000000))
  { 
    cat(paste0(gsub(pattern=base_dir,replacement = "",fmri.filelist[which(filesizes<10000000)]),"\n"))  
    prompt = utils::menu(c("Yes", "No"), title="\nThe above fMRI volumes have filesizes <10MB. Continue anyway?")
    if (prompt==2)
    {
      opt <- options(show.error.messages=FALSE)
      on.exit(options(opt))
      stop()
    }
  }
  
  ##subject level checks
  if(dir.check==T)
  {
    fmri.dir.list=unique(dirname(fmri.filelist))
    all_dir.check=rep(NA, length(dir.list))
    sub_dir.check=rep(NA, length(sub.list))
    fmri_dir.check = matrix(NA, nrow = length(fmri.dir.list), ncol = 2)
    
    #check for subject directories without at least one fMRI volume
    for (subj in 1:length(sub.list))  {sub_dir.check[subj]=length(which(stringr::str_detect(pattern = sub.list[subj],string = fmri.filelist)==T))}
    
    if(any(sub_dir.check==0)) 
    {
      cat(paste0("\nThe folowing subject directories do not contain at least 1 fMRI voluume; they will be ignored:\n"))
      cat(paste0(gsub(pattern=base_dir,replacement = "",sub.list[sub_dir.check==0]),"\n"))
      exc.sub=sub.list[sub_dir.check==0]
      sub.list=sub.list[-which(sub_dir.check==0)]
    }
    
    #check for fMRI directories without at least one fMRI volumes
    for (dir in 1:length(dir.list))  {all_dir.check[dir]=length(which(stringr::str_detect(pattern = dir.list[dir],string = fmri.filelist)))}
    
    if(any(all_dir.check==0)) 
    {
      dir.empty=dir.list[all_dir.check==0]
      
      if(any(sub_dir.check==0)) #remove invalid subjects detected earlier
      {
        remove.idx=list()
        for (sub.exc in 1:length(exc.sub))  {remove.idx[[sub.exc]]=which(stringr::str_detect(pattern = exc.sub[sub.exc],string = dir.empty)==T)}
        
        remove.idx=unlist(remove.idx)
        remove.idx=remove.idx[!is.na(remove.idx)]
        
        if(length(remove.idx>0))  {dir.empty=dir.empty[-remove.idx]}
      }
      
      if(length(dir.empty>0))
      {
        cat(paste0("\nThe folowing subjects' fMRI directories do not contain an fMRI volume and will be ignored:\n"))
        cat(paste0(gsub(pattern=base_dir,replacement = "",dir.empty),"\n"))     
      }
    }
    
    #check for missing movement files or fMRI directories containing more than one fMRI volume
    for (fmri_dir in 1:length(fmri.dir.list))
    {
      fmri_dir.check[fmri_dir,1]=length(which(stringr::str_detect(pattern = fmri.dir.list[fmri_dir],string = movement.filelist)==T))
      fmri_dir.check[fmri_dir,2]=length(which(stringr::str_detect(pattern = fmri.dir.list[fmri_dir],string = fmri.filelist)==T))
    }
    
    if(any(fmri_dir.check[,1]==0))
    {
      cat("\nThe following subject fMRI directories do not contain a movement file:\n")
      cat(paste0(gsub(pattern=base_dir,replacement = "",fmri.dir.list[fmri_dir.check[,1]==0]),"\n"))
      stop()
    }
    if(any(fmri_dir.check[,2]>1))
    {
      cat("\nThe following subject fMRI directories contain more than one fMRI volume:\n")
      cat(paste0(gsub(pattern=base_dir,replacement = "",fmri.dir.list[fmri_dir.check[,2]>1]),"\n"))
      stop()
    }
    cat("\n Directory structure check completed.\n\n")
  }
  
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
  
  ## loop thru subjects
  cat(paste("\nProcessing",length(sub.list), "valid subjects out of",N.orig,"listed subject directories...\n\n",sep=" "))
  for (sub in 1:NROW(sub.list))
  {
    cat(paste(sub.list[sub]," ",sub,"/",length(sub.list),"...",sep=""))
    if(((!file.exists(paste(output_dir,"/headmotion/",sub.list[sub],".txt",sep="")) | !file.exists(paste(output_dir,"/FCmat/",sub.list[sub],".csv",sep=""))) & overwrite==F) | overwrite==T)
    {
      start=Sys.time()
      #filepaths
      headmotion.param=list()
      headmotion.param$subject=sub.list[sub]
      fmri.path=fmri.filelist[which(stringr::str_detect(pattern = sub.list[sub],string = fmri.filelist)==T)]
      movement.path=movement.filelist[which(stringr::str_detect(pattern = sub.list[sub],string = movement.filelist)==T)]
      
      ##check number of movement files
      if(length(movement.path)>1)
      {
        movement.dat=list()
        frame.ends=list()
        frame.starts=list()
        
        for(volume in 1:length(movement.path))
        {
          mov.dat=read.table(movement.path[volume])
          if(NCOL(mov.dat)==1)  
          {movement.dat[[volume]]=mov.dat$V1} else if(NCOL(mov.dat)>=6)  
          {movement.dat[[volume]]=mov.dat[,1:6]} else  
          {movement.dat[[volume]]=NA}
          
          ##if its not the last volume 
          if(volume!=length(movement.path))  {frame.ends[[volume]]=c(NROW(movement.dat[[volume]]),NROW(movement.dat[[volume]])+1)}
          remove(mov.dat)
        }
        
        ##edit frame end indices only if there are more than 2 volumes; i.e. frame ends + number of frames in previous volumes
        if(length(movement.path)>2)  
        {
          for(volume in 2:(length(movement.path)-1))  
          {
            frame.starts[[volume]]=frame.ends[[volume-1]][2]+frame.ends[[volume]][1]
            frame.ends[[volume]]=frame.ends[[volume-1]][1]+frame.ends[[volume]]
          }
          frame.starts=c(frame.ends[[1]][2],unlist(frame.starts))
        }
        
        if(NCOL(movement.dat[[1]])==1)  {movement.dat=unlist(movement.dat)
        } else if (NCOL(movement.dat[[1]])>1)
        {
          movement.dat=do.call(rbind, movement.dat)
          movement.reg=movement.dat
          movement.dat=extractFD(movement.dat)
        }
        
        frame.ends=unlist(frame.ends)
        frame.start.dat=rep(0,NROW(movement.dat))
        frame.start.dat[frame.starts]=1
        
      } else
      {
        mov.dat=read.table(movement.path)
        if(NCOL(mov.dat)==1)  {movement.dat=mov.dat$V1
        } else if(NCOL(mov.dat)>=6)  
        {movement.dat=extractFD(mov.dat[,1:6])
        } else  
        {movement.dat=NA}
      }
      if(!any(is.na(movement.dat)))  {headmotion.param$mean=mean(movement.dat)
      } else  {headmotion.param$mean="Unable to compute. Missing/corrupted head motion data"}
      
      ##check number of fmri volumes, if multiple fmri volumes are detected, they are to be concatenated
      if(length(fmri.path)>1)
      {
        xii.list=list()
        for(vol in 1:length(fmri.path))  {xii.list[[vol]]=ciftiTools::read_xifti(fmri.path[vol], brainstructures="all")}
        xii.base=xii.list[[1]]
        
        for(vol in 1:length(fmri.path))
        {
          if(vol==1)  {xii.mat=scale(as.matrix(xii.list[[1]]))
          } else  {xii.mat=cbind(xii.mat,scale(as.matrix(xii.list[[vol]])))}
        }
        remove(xii.list,vol)
      } else
      {
        xii.base=ciftiTools::read_xifti(fmri.path, brainstructures="all")
        xii.mat=scale(as.matrix(xii.base))
      }
      headmotion.param$total_frames=NCOL(xii.mat)
      
      ##band pass filter
      if(length(bandpass)==2) { xii.mat=t(data.matrix(sapply(data.frame(t(xii.mat)), bandpassfilter, band=bandpass,TR=TR)))}
      
      ##scrubbing
      if((scrub==T & length(which(movement.dat>motion.thresh))!=0)) #if there are no frames less then the motion threshold, no scrubbing is needed
      {
        #identify frames for scrubbing
        exc_frames=which(movement.dat>motion.thresh)
        exc_frames=exc_frames[order(exc_frames)]
        exc_framesOLD=exc_frames
        
        if(length(exc_framesOLD)>1) #if there is only 1 bad frame, there is no between-frames interval to scrub
        {
          for (frameno in 1:(NROW(exc_framesOLD)-1))
          {
            if((exc_framesOLD[frameno+1]-exc_framesOLD[frameno]) <5)
            {
              if(length(movement.path)>1)
              {
                totest=match(frame.ends,(exc_framesOLD[frameno]:exc_framesOLD[frameno+1]))
                totest=totest[-is.na(totest)]
                if(length(totest)!=2)  {exc_frames=c(exc_frames,exc_framesOLD[frameno]:exc_framesOLD[frameno+1])}
              } else  {exc_frames=c(exc_frames,exc_framesOLD[frameno]:exc_framesOLD[frameno+1])}
            }
          }
          exc_frames=unique(exc_frames)
          exc_frames=exc_frames[order(exc_frames)]
        } 
        
        ##proceed with generating FC matrices if there are at least 10 frames
        if((headmotion.param$total_frames-length(exc_frames))>10)
        {
          #remove frames if necessary
          #if there is any missing movement data, or number of timepoints in the headmotion file do not match the number of frames, no scrubbing will be carried out.
          if(NCOL(xii.mat)!=length(movement.dat) | any(is.na(movement.dat)))
          {
            cat(paste("Missing timepoints in one of the movement files.\nScrubbing will not be carried out for",sub.list[sub]))
            xii.scrubbed=xii.mat
            movement.reg.scrubbed=movement.reg
            headmotion.param$no_frames_removed="Missing timepoints in one of the movement files., scrubbing cannot be carried out"
          } else
          {
            xii.scrubbed=xii.mat[,-exc_frames]
            movement.reg.scrubbed=movement.reg[-exc_frames,]
            frame.start.dat=frame.start.dat[-exc_frames]
            headmotion.param$no_frames_removed=length(exc_frames)
          }
          n_timepoints=NCOL(xii.scrubbed)
          
        ##confounds regression
          #selecting confounds for confounds regression
          if(motion.reg==T & GSR==F)  {confounds=movement.reg.scrubbed
          } else if(motion.reg==T & GSR==T) {confounds=cbind(movement.reg.scrubbed,colMeans(xii.scrubbed))
          } else if(motion.reg==F & GSR==T) {confounds=colMeans(xii.scrubbed)}
          
          if(GSR==T | motion.reg==T)
          {
            #computing derivatives and squares of confounds
            nconfound=NROW(confounds)
            confounds=data.matrix(confounds)
            confounds.derivative=confounds
            confounds.derivative[2:nconfound,]=confounds[2:nconfound,]-confounds[1:(nconfound-1),]
            confounds.derivative[1,]=0
            if(exists("frame.start.dat")){confounds.derivative[which(frame.start.dat==1),]=0}
            all.confounds=cbind(confounds, confounds^2,confounds.derivative,confounds.derivative^2)
            
            #regressing out confounds
            dmat=xii.scrubbed
            dmat=fMRItools::nuisance_regression(dmat, cbind(1, all.confounds))
            xii.regressed=ciftiTools::newdata_xifti(xii.base, dmat)
            xii.final=ciftiTools::move_from_mwall(xii.regressed, NA)
            
            remove(dmat,xii.base,xii.scrubbed,xii.regressed,confounds,confounds.derivative,nconfound)
          } else
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
          #reorder parcel indices for visualization purpose
          xii_pmean=xii_pmean[,c(1:atlas,reorder.subcortical.idx)]
          
          #output headmotion stats
          headmotion.row=matrix(c(headmotion.param$subject,headmotion.param$mean,headmotion.param$total_frames,headmotion.param$no_frames_removed),nrow=1)
          write.table(headmotion.row, file=paste(output_dir,"/headmotion/",sub.list[sub],".txt",sep=""), row.names = F, col.names = F, sep=",",quote=F)

          #output FC matrix/timeseries
          if (output.timeseries==F)  {write.table(cor(xii_pmean), file=paste(output_dir,"/FCmat/",sub.list[sub],".csv",sep=""), row.names = F, col.names = F, sep=",",quote=F)
          }  else if(output.timeseries==T) {write.table(xii_pmean, file=paste(output_dir,"/FCmat/",sub.list[sub],".csv",sep=""), row.names = F, col.names = F, sep=",",quote=F)}
          
          remove(xii.final, sub_keys,brain_vec,xii_pmean,data_p, headmotion.row, headmotion.param)
          end=Sys.time()
          
          cat(paste(" completed in",round(difftime(end,start, units="secs"),1),"secs\n",sep=" "))   
        } else  {cat("Less than 10 valid frames exist after scrubbing. FC matrix will not be generated for this subject.\n")}
      }
    }  else  {cat(" Not processed; output files from a previous run are detected.\n")}
  }
  
  ##compiling headmotion parameters across subjects
  cat("\n Compiling report...")  
  headmotion.output.filelist=list.files(paste0(output_dir,"/headmotion/"))
  if(scrub==T)  
  {
    report.all=data.frame(matrix(NA,nrow=NROW(headmotion.output.filelist),ncol=4))
    colnames(report.all)=c("subject","mean FD", "total frames","no.frames_removed")
  } else
  {
    report.all=data.frame(matrix(NA,nrow=NROW(headmotion.output.filelist),ncol=3))
    colnames(report.all)=c("subject","mean FD", "total frames")
  }
  
  for (sub in 1:length(headmotion.output.filelist))
  {
    report.all[sub,]=read.table(paste0(output_dir,"/headmotion/",headmotion.output.filelist[sub]),header = F,sep=",")[1,]
  }
  write.table(report.all,file=paste0(output_dir,"/report.csv"),row.names = F,sep=",")
  cat(paste0("done. Report written to ",output_dir,"/report.csv")) 
}

########################################################################################################
########################################################################################################
##function to extract FD, adapted from fMRItools::FD

extractFD=function(mot_dat)
{
  detrend = F
  mot_dat[, 4:6] <- mot_dat[, 4:6]
  if (!isFALSE(detrend)) {
    if (isTRUE(detrend)) {
      detrend <- 4
    }
    mot_dat <- fMRItools::nuisance_regression(mot_dat, cbind(1, fMRItools::dct_bases(nrow(mot_dat), detrend)))
  }
  mot_datdiff <- apply(mot_dat, 2, diff, lag = 1)
  return(c(rep(0, 1), rowSums(abs(mot_datdiff))))
}
########################################################################################################
########################################################################################################
##band pass filter function

bandpassfilter=function(data,TR,band)
{
  nframes = length(data)
  t = seq(0, nframes, by = TR)
  
  fft_data <- fft(data)
  freq <- (0:(nframes-1)) / (nframes * TR)
  
  fft_data[abs(freq) > band[2]] <- 0
  fft_data[abs(freq) < band[1]] <- 0
  
  return(Re(fft(fft_data, inverse = TRUE)) / nframes)
}
