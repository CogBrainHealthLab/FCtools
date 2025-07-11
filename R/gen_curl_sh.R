########################################################################################################
########################################################################################################
#' @title gen_curl_sh
#'
#' @description function to generate curl shellscript for downloading data from openneuro
#'
#' @details This function searches the `sh` input file for keywords matching strings specified in `include_or`, `include_and`, `exclude` and `subjects`, and returns a vector of curl commands as well as output a .sh file
#'
#' @param sh The filepath to sh file containing the curls commands as generated by openneuro
#' @param include_and The curl commands have to contain all of the keywords specfied here will be included
#' @param include_or The curl commands have to contain at least one of the keywords specfied here will be included
#' @param exclude The curl commands which contain at least one of the keywords specfied here will be excluded
#' @param subjects If specified only curl commands that include these subjects will be included
#' @param filename filename of the .sh output file. Set to `curl_cmd.sh` by default.
#' @returns a vector of curl commands
#'
#' @importFrom stringr str_detect
#' @export
########################################################################################################
########################################################################################################
gen_curl_sh=function(sh,include_and,include_or,exclude, subjects,filename="curl_cmd.sh")
{
  curl_cmd=read.table(sh,sep = "\n")$V1
  dataset=curl_cmd[which(stringr::str_detect(string=curl_cmd,pattern="dataset_description.json"))]
  curl_cmd=curl_cmd[-which(stringr::str_detect(string=curl_cmd,pattern="dataset_description.json"))]
  if(!missing(include_and))
  {
    for(keyword in include_and)
    {
    curl_cmd=curl_cmd[which(stringr::str_detect(string=curl_cmd,pattern=keyword))]
    }
  }
  if(!missing(include_or))
  {
    for(keyword in 1:length(include_or))
    {
      if(keyword==1)
      {
      curl_cmd=curl_cmd[which(stringr::str_detect(string=curl_cmd,pattern=include_or[keyword]))]
      } else
      {
      curl_cmd=c(curl_cmd,curl_cmd[which(stringr::str_detect(string=curl_cmd,pattern=include_or[keyword]))])
      }
    }
  }
  if(!missing(exclude))
  {
    for(keyword in exclude)
    {
    curl_cmd=curl_cmd[-which(stringr::str_detect(string=curl_cmd,pattern=keyword))]
    }
  }
  if(!missing(subjects))
  {
    for(sub in 1:length(subjects))
    {
      if(sub==1)
      {
      curl_cmd=curl_cmd[stringr::str_detect(string=curl_cmd,pattern =subjects[sub])]
      } else
      {
      curl_cmd=c(curl_cmd,curl_cmd[stringr::str_detect(string=curl_cmd,pattern =subjects[sub])])
      }
    }
  }
  write.table(c(dataset,curl_cmd),file=filename,quote = F, row.names = F, col.names = F)
  return(c(dataset,curl_cmd))
}
