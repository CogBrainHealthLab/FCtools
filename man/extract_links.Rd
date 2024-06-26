% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCP_dat_func.R
\name{extract_links}
\alias{extract_links}
\title{extract_links}
\usage{
extract_links(
  manifest = "datastructure_manifest.txt",
  files,
  filename = "downloadlist.txt"
)
}
\arguments{
\item{manifest}{the filepath of the manifest file. Set to \code{datastructure_manifest.txt} by default.}

\item{files}{a vector of keywords that contained within the file names of the files to be downloaded.}

\item{filename}{The desired filename, with a *.txt extension, of the text file containing the links. Set to \code{download_list.txt} by default}
}
\value{
outputs a .txt file containing the links filtered from the data manifest file
}
\description{
Extracting links from data manifests
}
\details{
This function extracts links from NDA formatted data manifest files
}
\examples{
\dontrun{
extract_links(manifest="datastructure_manifest.txt", filename="download_list.txt", files=c("REST_Atlas_MSMAll_hp0_clean.dtseries.nii", "Movement_RelativeRMS.txt"))
}
}
