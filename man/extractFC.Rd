% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCP_dat_func.R
\name{extractFC}
\alias{extractFC}
\title{extractFC}
\usage{
extractFC(
  wb_path,
  base_dir = "fmriresults01/",
  atlas,
  subjects,
  motion.thresh = 0.25,
  scrub = F,
  GSR = F,
  task,
  extension = "_Atlas_MSMAll_hp0_clean.dtseries.nii",
  movement.extension = "Movement_RelativeRMS.txt",
  output_dir = "output_dir",
  overwrite = TRUE,
  dir.check = TRUE
)
}
\arguments{
\item{wb_path}{The filepath to the workbench folder that you have previously downloaded and unzipped}

\item{base_dir}{The filepath to directory containing the subject folders. Set to \verb{fmriresults01/"} by default.}

\item{atlas}{The version (number of parcels) of the Schaefer atlas; it has to be in multiples of 100, from 100 to 1000. The specified atlas template will be automatically downloaded from here if it does not exist in the current directory.}

\item{motion.thresh}{the threshold used for scrubbing frames. Set to 0.25 by default.}

\item{scrub}{If set to TRUE, remove frames with excessive head motion (relative RMS displacement> 0.25) and also remove segments in between two time-points if the segment has less than 5 frames. This is the ‘full scrubbing; method described in \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811913009117)}{Power et. al (2014)}. This is not recommended for task-based fMRI volumes.}

\item{GSR}{If set to TRUE, global signal and its square-term will be regressed out from the fMRI timeseries data.}

\item{task}{The name of the task (case-sensitive), without any numbers or acquisition direction labels}

\item{extension}{The filename extension of the fMRI volumes. Set to \verb{_Atlas_MSMAll_hp0_clean.dtseries.nii} by default}

\item{movement.extension}{The filename extension of the head motion files. Set to \code{Movement_RelativeRMS.txt} by default. Note the head motion file should be a text file containing either a single column of numbers (RMS measurements), or at least 6 columns (3 displacement + 3 rotational vectors). If there are more than 6 columns, only the first 6 columns will be used to compute framewise displacement.}

\item{output_dir}{The output directory where the FC matrices and headmotion parameters will be saved. A default \code{output_dir} directory will be created if it does not exist.}

\item{overwrite}{If set to \code{FALSE}, subjects that already have their FC matrices in \code{output_dir} will be skipped. Set to \code{TRUE} by default.}

\item{dir.check}{check the directory structure to ensure that each fMRI directory should only contain one fMRI file and at least one movement file. Set to \code{TRUE} by default.}
}
\value{
outputs M x M matrices in \verb{output_dir\\FCmat}, headmotion measurements in \verb{output_dir\\headmotion} and a summary headmotion measurements across all subjects \verb{output_dir\\report.csv}.
}
\description{
Extracting FC matrices from HCP-preprocessed fMRI volumes
}
\details{
This function extracts FC matrices from fMRI volumes preprocessed with the HCP pipeline with the option
of carrying out global signal regression and scrubbing. The directory structure and file names must adhere to the NDA/HCP data structure
and file naming conventions. \href{https://www.humanconnectome.org/software/get-connectome-workbench)}{Connectome Workbench} has to be
installed prior to running the function.
}
\examples{
\dontrun{
extractFC(wb_path = "/home/junhong.yu/workbench/",task="REST",atlas=100,GSR=TRUE,scrub=TRUE,output_dir = "FCmat")
}
}
