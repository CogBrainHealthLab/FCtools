############################################################################################################################
############################################################################################################################
#' @title extract_headmotion
#'
#' @description Extracting head motion measurements.
#'
#' @details This function extracts head motion measurements from an fMRIprep output directory, and outputs these measurements in a .csv file
#'
#' @param filename the desired filename, with a *.csv extension, of the output.Set to 'motiondat.csv' by default
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
extract_headmotion=function(filename="motiondat.csv")
{
  filelist=list.files(pattern=".html")
  subjlist=gsub(".html","",filelist)
  motiondat=data.frame(cbind(subjlist,rep(0,NROW(subjlist)),rep(0,NROW(subjlist)),rep(0,NROW(subjlist)),rep(0,NROW(subjlist)),rep(0,NROW(subjlist))))
  colnames(motiondat)=c("subj","no.frames","FD","RMSD","FD.20","RMSD.25")
  filename=gsub(subjlist[1],"",list.files(path=paste(subjlist[1],"/func/", sep=""),pattern="confounds_timeseries.tsv"))
  for (subjno in 1:NROW(subjlist))
  {
    confounds=read.table(file = paste(subjlist[subjno],"/func/",subjlist[subjno],filename,sep=""), sep = '\t', header = TRUE)
    confounds[confounds=="n/a"]=NA

    motiondat$RMSD[subjno]=mean(as.numeric(confounds$rmsd), na.rm=T)
    motiondat$no.frames=NROW(confounds)
    motiondat$FD[subjno]=mean(as.numeric(confounds$framewise_displacement), na.rm=T)
    motiondat$RMSD.25[subjno]=NROW(which(confounds$rmsd>0.25))/nrow(confounds)
    motiondat$FD.20[subjno]=NROW(which(confounds$framewise_displacement>0.20))/nrow(confounds)
  }
  write.csv(motiondat,file=filename)
}
########################################################################################################
########################################################################################################

#' @title vectorizeFC
#'
#' @description vectorizing and concatenating and M x M matrices
#'
#' @details This function reads the M x M matrices (saved as .csv files) in a directory, vectorizes and concatenates them into a single N x ((M x M)-M)/2 matrix. N= number of subjects; M=number of nodes in the connectome
#'
#' @param FC_dir the directory containing the M x M matrices (saved as .csv files)
#' @param filename the output format of the concatenated group level matrix. Only filenames with the *.csv or *.rds extensions are accepted. The .rds format uses less disk space. Set to 'FCmat.rds' by default
#' @param fisherz `TRUE` or `FALSE` option specifying if the edge values should undergo a fisher-Z transformation. This would be applicable in cases and functional connectivity is calculated as a pearson's correlation coefficient. Set to `TRUE` by default.
#' @param timeseries `TRUE` or `FALSE` option specifying if the parcellated timeseries are used as inputs. Set to `FALSE` by default.
#' @returns Depending on whether a filename with a .csv or .rds file extension is specified, outputs a .csv file containing the  N x ((M x M)-M)/2 matrix or a .rds file containing a list object: 1)N x ((M x M)-M)/2 matrix, 2) a vector of subject IDs.
#'
#' @examples
#' \dontrun{
#' vectorizeFC(FC_dir="FCmat", filename="FCmat.rds", fisherz=TRUE)
#' }
#' @importFrom stringr str_detect
#' @importFrom psych fisherz
#' @export
#'
#'
########################################################################################################
########################################################################################################
vectorizeFC=function(FC_dir, filename="FCmat.rds", fisherz=TRUE, timeseries=FALSE)
{
  ##check filename extension
  if(stringr::str_detect(string=filename, pattern=".csv")==F & stringr::str_detect(string=filename, pattern=".rds")==F)
  {
    stop("invalid filename extension")
  }

  subjlist=list.files(path=FC_dir)

  if(timeseries==F)
  {
    FCformat=read.csv(paste(FC_dir,"/",subjlist[1], sep=""), header=F)

    #if subjects' FCmat is already vectorized
    if(dim(FCformat)[1]==1)
    {
      FCmat=matrix(NA,nrow=NROW(subjlist), ncol=dim(FCformat)[2])
      for(sub in 1:NROW(subjlist))
      {
        FCmat[sub,]=unlist(read.table(paste(FC_dir,"/",subjlist[sub], sep=""),header = F,sep=",",col.names = F))
        cat(paste(subjlist[sub],"\n"))
      }
      if(fisherz==T)
      {
        FCmat=psych::fisherz(FCmat)
      }
    }

    #if subjects' FCmat has not been vectorized
    if(dim(FCformat)[1]>1)
    {
      x=dim(FCformat)[1]
      FCmat=matrix(NA,nrow=NROW(subjlist), ncol=((x*x-x)/2))
      for(sub in 1:NROW(subjlist))
      {
        FCmat.temp=read.table(paste(FC_dir,"/",subjlist[sub], sep=""),header = F,sep=",")
        FCmat[sub,]=FCmat.temp[upper.tri(FCmat.temp, diag = FALSE)]
        cat(paste(subjlist[sub],"\n"))
        remove(FCmat.temp)
      }
      if(fisherz==T)
      {
        FCmat=psych::fisherz(FCmat)
      }
    }
  } else if(timeseries==T)
  {
    FCformat=read.csv(paste(FC_dir,"/",subjlist[1], sep=""), header=F)
    x=NCOL(FCformat)
    FCmat=matrix(NA,nrow=NROW(subjlist), ncol=((x*x-x)/2))
    for(sub in 1:NROW(subjlist))
    {
      FCmat.temp=read.table(paste(FC_dir,"/",subjlist[sub], sep=""),header = F,sep=",")
      FCmat[sub,]=FCmat.temp[upper.tri(cor(FCmat.temp), diag = FALSE)]
      cat(paste(subjlist[sub],"\n"))
      remove(FCmat.temp)
    }
    if(fisherz==T)
    {
      FCmat=psych::fisherz(FCmat)
    }
  }

  ##output file
  if(stringr::str_detect(string=filename, pattern=".csv"))
  {
    write.table(FCmat,file = filename,row.names = F,col.names = F, sep=",")
    cat(paste("Concatenated FC matrices saved to",filename,"\n",sep=" "))
  } else if(stringr::str_detect(string=filename, pattern=".rds"))
  {
    saveRDS(list(FCmat,gsub(".csv","",subjlist)),file = filename)
    cat(paste("Concatenated FC matrices saved to",filename,"\n",sep=" "))
  }
}

########################################################################################################
########################################################################################################
