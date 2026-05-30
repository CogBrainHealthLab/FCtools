############################################################################################################################
############################################################################################################################
#' @title edges_to_networks
#'
#' @description convert an edge level FC matrix into a network level FC matrix
#'
#' @details This function first identifies the unique network pairs in the appropriate FC atlas and then averages the edges within each of the network pairs
#' @param FCmat an FC matrix or vector
#' 
#' @returns returns a network level FC matrix
#'
#' @examples
#' edges_to_networks(runif(23871))
#' 
#' @export
############################################################################################################################
############################################################################################################################
edges_to_networks=function(FCmat)
{
  if(!any(class(FCmat)=="matrix"))  {FCmat=t(data.matrix(FCmat))}
  
  FCmat.n_nodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FCmat))) / (2 * 0.5)
  atlas.no=match(FCmat.n_nodes,c(90,119,219,246))
  if(is.na(atlas.no)) {stop(paste0("FCmat contains ", NCOL(FCmat),"columns, which is consistent with the support atlases"))}
  
  labels=get("labels_dat")[[atlas.no]]
  
  ##generate unique network pairs and their labels
  network.pairs=as.matrix(expand.grid(1:max(labels$region), 1:max(labels$region)))
  network.pairs=network.pairs[network.pairs[,1] <= network.pairs[,2], ]
  
  network_to_no=unique(cbind(labels$regionlabel,labels$region))
  network.pairs.labs=network.pairs
  
  for(network.no in as.numeric(network_to_no[,2]))
  {
    network.pairs.labs=matrix(gsub(network.no,network_to_no[match(network.no,as.numeric(network_to_no[,2])),1],network.pairs.labs),ncol=2)
  }
  
  netFC=matrix(NA, nrow = NROW(FCmat),ncol=36)
  FCmat.template=matrix(NA,nrow=219, ncol=219)
  FCmat.template[upper.tri(FCmat.template)]=1:23871
  for(edge in 1:NROW(network.pairs))
  {
    ##identifying edge indices belong to a network
    edge.idx=unique(as.numeric(FCmat.template[which(labels$region==network.pairs[edge,1]),which(labels$region==network.pairs[edge,2])]))
    edge.idx=edge.idx[!is.na(edge.idx)]
    
    ##average edges across network groupings
    if(NROW(FCmat)>1) {netFC[,edge]=rowMeans(FCmat[,edge.idx])}
    else {netFC[,edge]=mean(FCmat[,edge.idx])}
  }
  colnames(netFC)=paste(network.pairs.labs[,2],"to", network.pairs.labs[,1])
  return(netFC)
}