## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title reorder.subcortical
#'
#' @description to reorder subcortical nodes for visualization
#'
#' @details The subcortical node elements in the previously extracted HCP matrices are not ordered correctly. Hence the visualizations (`vizHeatmap()` or `vizConnectogram()`) of the subcortical nodes generated from these HCP matrices are not accurate; Specifically the positioning of the subcortical nodes within the subcortical network is wrong, while the cortical nodes are unaffected. This function thus corrects the ordering of the subcortical nodes
#'
#' @param data a vector of edge values with a length of 7021 or 23871
#'
#' @returns a vector of edge values with the correct subcortical nodes ordering
#'
#' @examples
#'
#' FC.vector=runif(23871, min = -1, max = 1)
#' reorder.subcortical(FC.vector)
#'
#' @export
##Main function
########################################################################################################
########################################################################################################
reorder.subcortical=function(data)
{
data=data.matrix(data)  
if (NCOL(data)==23871)  {n_nodes=219
} else if (NCOL(data)==7021)  {n_nodes=119}

cort.nodes=n_nodes-19
reorder.subcortical.idx=c(1:cort.nodes,c(9,18,8,17,6,3,13,1,11,10,19,7,16,5,15,4,14,2,12)+cort.nodes)
FCmat.temp=matrix(0, nrow=n_nodes, ncol=n_nodes)
FCmat.temp[upper.tri(FCmat.temp,diag = F)]=data
FCmat.temp=FCmat.temp+t(FCmat.temp)
FCmat.reordered=FCmat.temp[reorder.subcortical.idx,reorder.subcortical.idx]
FC.vector.reordered=FCmat.reordered[upper.tri(FCmat.reordered,diag = F)]

return(FC.vector.reordered)
}
