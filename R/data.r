#' Labels data
#'
#' A list containing four data frames, listing atlas labels, index numbers, hemispheres and lobes they correspond to, for each available atlas (1=AAL-90 atlas, 2=Schaefer-100 atlas + 19 subcortical regions, 3=Schaefer-200 atlas + 19 subcortical regions, 4=Brainnetome atlas). 
#'
#' @name labels_dat
#' @docType data 
#' @keywords internal 
#' @format ## `labels_dat`
#' A list of four data.frames: () 
#' \describe{
#'   \item{vertices}{data frame with 90 rows (atlas parcellations), 6 columns (atlas label, region number, index order (old), index order (new), hemisphere (1=L,2=R), region label (anatomical area name) )}
#'   \item{vertices}{data frame with 119 rows (atlas parcellations), 6 columns (atlas label, region number, index order (old), index order (new), hemisphere (1=L,2=R), region label (anatomical area name) )}
#'   \item{vertices}{data frame with 219 rows (atlas parcellations), 6 columns (atlas label, region number, index order (old), index order (new), hemisphere (1=L,2=R), region label (anatomical area name) )}
#'   \item{vertices}{data frame with 246 rows (atlas parcellations), 6 columns (atlas label, region number, index order (old), index order (new), hemisphere (1=L,2=R), region label (anatomical area name) )}
#' }
'labels_dat'

#' Demo data 
#' 
#' A small sample of 4 connectivity vector derived from 219x219 FC matrices generated using the Schaefer-200 atlas + 19 subcortical regions from the freesurfer subcortical segmentations. 
#' 
#' @name demomat
#' @docType data 
#' @keywords internal 
#' @format ## `demomat`
#' \describe{matrix with 4 rows (number of individuals), 23871 columns (edge weights)}
'demomat'

#' Fsaverage5 
#' 
#' A surface template in fsaverage5 space for the glass brain plot
#' 
#' @name fs5brain
#' @docType data 
#' @keywords internal 
#' @format ## `fs5brain`
#' A list of two matrices: () 
#' \describe{
#'   \item{vertices}{matrix with 3 rows (XYZ coordinates), 20484 columns (number of vertices)}
#'   \item{edges}{matrix with 3 rows (triangle neighbours), 40960 columns (unique triangles)}
#'   }
'fs5brain'