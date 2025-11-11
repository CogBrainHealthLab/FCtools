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
  #read manifest and remove first row; the first row contains description of the column, hence not used.
  manifest=read.table(manifest,header = T)[-1,]

  #select subjects
  if(!missing(subjects))
  {
    sub.list=list()
    for(sub in 1:length(subjects))
    {
      sub.list[[sub]]=which(stringr::str_detect(string=manifest$associated_file,pattern =subjects[sub]))
    }
  filelist.incsub=manifest$associated_file[unique(unlist(sub.list))]
  
  } else
  {
    filelist.incsub=manifest$associated_file
  }
  
  if(!missing(task))
  {
    #select required task-related files

    files=c(paste0("_",task,"_.._Atlas_MSMAll.dtseries.nii"),
            paste0("_",task,"._.._Atlas_MSMAll.dtseries.nii"),
            paste0("_",task,".._.._Atlas_MSMAll.dtseries.nii"),
            "brainmask_fs.2.nii.gz",
            "_SBRef.nii.gz",
            "Movement_Regressors.txt",
            "Movement_AbsoluteRMS.txt",
            paste0("_",task,".._...nii.gz"),
            paste0("_",task,"._...nii.gz"),
            paste0("_",task,"_...nii.gz"))
    
    idx.list=list()
    for(file in 1:NROW(files))
    {
      idx.list[[file]]=which(stringr::str_detect(pattern = files[file],string = filelist.incsub)==T)
    }
    
    filelist.alltasks=filelist.incsub[unique(unlist(idx.list))]

    #filter tasks
    task.idx.list=list()
    for(t in 1:length(task))
    {
      task.idx.list[[t]]=which(stringr::str_detect(pattern = task[t],string = filelist.alltasks)==T)
    }
    filelist.sel.tasks=filelist.alltasks[unique(unlist(task.idx.list))]
    
    #add MNINonlinear files
    MNIfiles=c("MNINonLinear/T1w.nii.gz",
               "MNINonLinear/ribbon.nii.gz",
               "MNINonLinear/aparc\\+aseg",
               "MNINonLinear/T1w.nii.gz",
               "MNINonLinear/brainmask_fs.nii.gz",
               "MNINonLinear/brainmask_fs.2.nii.gz")
    MNI.idx.list=list()
    for(file in 1:length(MNIfiles))
    {
      MNI.idx.list[[file]]=which(stringr::str_detect(pattern = MNIfiles[file],string = filelist.incsub)==T)
    }
    filelist.MNI=filelist.incsub[unique(unlist(MNI.idx.list))]
    filelist.sel.tasks.MNI=c(filelist.sel.tasks,filelist.MNI)
  }
  #add surface files
  if(!missing(surf))
  {
    surf_files=c("lh.cortex.label",
            "rh.cortex.label",
            "lh.sphere.reg",
            "rh.sphere.reg",
            paste0("lh.",surf),
            paste0("rh.",surf),
            "aseg.stats")
    surf.idx.list=list()
    for(file in 1:NROW(surf_files))
    {
      surf.idx.list[[file]]=which(stringr::str_ends(pattern = surf_files[file],string = filelist.incsub)==T)
    }
    filelist.surf=filelist.incsub[unique(unlist(surf.idx.list))]

    if(!missing(task))
      {
      filelist.sel.tasks.MNI=c(filelist.surf,filelist.sel.tasks.MNI)
      }
    else
      {
      filelist.sel.tasks.MNI=filelist.surf
      }
  }
  #check if files were found
  if(length(filelist.sel.tasks.MNI)==0)
  {
    stop("No files were found, check your task name")
  }
  #output filelist as a text file
  write.table(filelist.sel.tasks.MNI[order(filelist.sel.tasks.MNI)],file=filename, quote = F, row.names = F, col.names = F)
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
