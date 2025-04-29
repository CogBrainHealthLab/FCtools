############################################################################################################################
############################################################################################################################
#' @title editJSON_IntendedFor
#'
#' @description A tool to edit .JSON files' IntendedFor entry so that the fmap files can be correctly identified and processed by fMRIprep. This currently does not work with BIDS directories that contain multiple sessions of data
#'
#' @details This function searches for .JSON files at the root of fMRIprep BIDS directory, edits their IntendedFor entry before overwriting the original JSON file
#'
#' @param json_files filename suffix for one (e.g., `"_phasediff.json"`) or multiple JSON files (e.g., `c("_phasediff.json", "magnitude1.json", "magnitude2.json")`)
#' @param nii_files one or more filename suffix of the nii or nii.gz files (e.g., `"_task-rest_bold.nii.gz"`) to be added to the IntendedFor entry of the specified JSON files. 
#' @param subjects If specified, will only processed the selected subjects
#' @param base_dir base directory of the *bold.nii.gz files. set to "func/" by default.
#'
#' @examples
#' \dontrun{
#' editJSON_IntendedFor(json_files="_phasediff.json",nii_files="_task-rest_bold.nii.gz")
#' }
#' @importFrom stringr str_detect
#' @importFrom rjson fromJSON toJSON
#' @export
########################################################################################################
########################################################################################################

editJSON_IntendedFor=function(json_files,nii_files,subjects,base_dir="func/")
{
  #if only a single json file is specified
  cat("\nSearching for .json files...")
  if(length(json_files)==1) 
  {
    filelist=list.files(recursive = T,pattern = json_files)
    if(length(filelist)==0) {stop(paste0("No *",json_files[json_file]," files are found"))}
    
    sublist=gsub(pattern =json_files,"",basename(filelist))
    sub.list.idx=list()
    
    #if subjects is specified, filter out excluded subjects
    if(!missing("subjects"))
    {
      for (subj in 1:length(subjects))  {sub.list.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = filelist)==T)}
      
      sub.list.idx=unlist(sub.list.idx)
      filelist=filelist[sub.list.idx]
      sublist=sublist[sub.list.idx]
    }
    cat("\nEditing .json files...")
    #loop across subjets to edit JSON files
    for (sub in 1:length(sublist))
    {
      json.obj=rjson::fromJSON(file=filelist[sub])
      json.obj$IntendedFor=paste0(base_dir,sublist[sub],nii_files)
      
      write(rjson::toJSON(json.obj,2),file=filelist[sub])
    }  
  } 
  
  #if multiple json files are specified
  if(length(json_files)>1) 
  {
    for (json_file in 1:length(json_files))
    {
    
      filelist=list.files(recursive = T,pattern = json_files[json_file])
      if(length(filelist)==0) {stop(paste0("No *",json_files[json_file]," files are found"))}
      
      sublist=gsub(pattern =json_files[json_file],"",basename(filelist))
      sub.list.idx=list()
      
      #if subjects is specified, filter out excluded subjects
      if(!missing("subjects"))
      {
        for (subj in 1:length(subjects))  {sub.list.idx[[subj]]=which(stringr::str_detect(pattern = subjects[subj],string = filelist)==T)}
        
        sub.list.idx=unlist(sub.list.idx)
        filelist=filelist[sub.list.idx]
        sublist=sublist[sub.list.idx]
      }
      
      #loop across subjets to edit JSON files
      for (sub in 1:length(sublist))
      {
        json.obj=rjson::fromJSON(file=filelist[sub])
        json.obj$IntendedFor=paste0(base_dir,sublist[sub],nii_files)
        
        write(rjson::toJSON(json.obj,2),file=filelist[sub])
      }  
    } 
  }
    cat("\nDone")
}
