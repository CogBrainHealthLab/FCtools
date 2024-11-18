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