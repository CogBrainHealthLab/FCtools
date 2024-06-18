## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title vizConnectogram
#'
#' @description Visualizing brain connectivity profiles with a connectogram
#'
#' @details This function takes a vector of edges and visualizes the edge-to-edge connectivity in a connectogram.
#'
#' @param data a vector of edge values with a length of 4005, 7021, 23871 or 30135.
#' @param hot color for the positive connections. Set to `#F8766D` by default.
#' @param cold color for the negative connections. Set to `#00BFC4` by default.
#' @param colorscheme an optional vector of color names or color codes to color code the networks.
#' @param filename output filename with a *.png file extension. Set to "conn.png" by default
#' @param edgethickness a value to adjust the thickness of the edges. Set to 0.8 by default.
#' @param colorbar `TRUE` or `FALSE` option on whether a colorbar for the connectivity strength will be displayed. Set to `TRUE` by default.
#' @returns outputs a .png image
#'
#' @examples
#'
#' results=sample(c(1,0, -1), 7021, replace = TRUE, prob = c(0.01, 0.98,0.01))
#' vizConnectogram(data=results, filename="FC_119.png")
#'
#' @importFrom ggplot2 ggplot aes scale_colour_gradient2 geom_point guides theme guide_colorbar element_text scale_color_manual expand_limits unit guide_legend
#' @importFrom igraph graph_from_adjacency_matrix get.edgelist E edge_attr
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar
#' @importFrom ggraph ggraph geom_edge_arc scale_edge_alpha_continuous scale_edge_color_manual geom_node_point geom_node_text theme_graph
#' @export

########################################################################################################
########################################################################################################


vizConnectogram=function(data, hot="#F8766D", cold="#00BFC4", edgethickness=0.8,filename="conn.png", colorscheme, colorbar=T)
{
  ##selecting atlas

  edge_lengths=c(4005,7021,23871,30135)

  if(is.na(match(length(data),edge_lengths)))
  {
    stop("The length of the input vector does not fit any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")
  } else
  {
    atlas=match(length(data),edge_lengths)
  }

  ##parameters for all atlases
  param=list(NA,NA,NA,NA,NA)
  names(param)=c("nodecol","nodesize","xlim","ylim","nodelevels")
  param$nodecol=list(c("#D53E4F","#FC8D59","#FEE08B","#FFFFBF","#E6F598","#99D594","#3288BD"),
                     c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD"))
  param$nodecol[3:4]=param$nodecol[2]
  param$nodesize=c(1.4,1.3,1,1)
  param$xlim=list(c(-1.15, 1.6),
                  c(-1.15, 1.5),
                  c(-1.2, 1.6),
                  c(-1.15, 1.5))
  param$ylim=list(c(-1.2, 1.2),
                  c(-1.2, 1),
                  c(-1.15, 1.2),
                  c(-1.2, 1))

  label=labels_dat[[match(length(data),edge_lengths)]]
  param$nodelevels=unique(label$regionlabel)
  label=label[order(label$oldorder),]
  label.neworder=label[order(label$neworder),]
  param$nodelevels=unique(label.neworder$regionlabel)
  label$regionlabel = factor(label$regionlabel,levels = param$nodelevels)
  if(missing("colorscheme")){colorscheme = param$nodecol[[atlas]]}
  nnodes=nrow(label)

  ##rehaping data into connectivity matrix

  conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
  conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = data
  conmat_NxN=conmat_NxNhalf+t(conmat_NxNhalf)

  nodeorder=as.numeric(rep(NA,nnodes))
  for (rowno in 1:nnodes){
    nodeorder[rowno]=which(label$neworder==rowno)
  }

  colnames(conmat_NxN)=label$labels
  rownames(conmat_NxN)=label$labels
  reordered=conmat_NxN[nodeorder,nodeorder]
  RegionsFC=as.factor(label$regionlabel)[nodeorder]

  ##graph object

  graphobjFC=igraph::graph_from_adjacency_matrix(reordered, mode="undirected", diag=FALSE, weighted=T)
  EvalFC=igraph::edge_attr(graphobjFC, "weight", index = igraph::E(graphobjFC))
  posnegFC=EvalFC
  posnegFC=replace(posnegFC, which(posnegFC < 0), "2_neg")
  posnegFC=replace(posnegFC, which(posnegFC!="2_neg"), "1_pos")

  igraph::edge_attr(graphobjFC, "weight", index = igraph::E(graphobjFC))=abs(EvalFC)
  igraph::edge_attr(graphobjFC, "posFC", index = igraph::E(graphobjFC))=posnegFC

  #plot

  FCplot=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
    ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness, show.legend = T) +
    ggraph::scale_edge_alpha_continuous(guide="none")+
    ggraph::scale_edge_color_manual(name="Edges", labels=c("Positive","Negative"),values=c(hot,cold))+
    ggplot2::scale_color_manual(values =colorscheme, name="Network")+
    ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=param$nodesize[atlas], shape=19,show.legend = T) +
    ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                        angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                        hjust = ifelse(x > 0, 0 ,1)), size=param$nodesize[atlas]) +
    ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
    ggplot2::expand_limits(x = param$xlim[[atlas]], y = param$ylim[[atlas]])+
    ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4),
                   legend.title=ggplot2::element_text(size=5,face = "bold"),
                   legend.text=ggplot2::element_text(size=5),
                   legend.position = c(1,0),legend.justification=c(1,0), legend.key.height = ggplot2::unit(c(0, 0, 0, 0), "cm"))

  if(colorbar==T)
  {
    FCplot=FCplot+
      ggplot2::guides(edge_color = ggplot2::guide_legend(override.aes = list(shape = NA)),color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))
  } else if(colorbar==F)
  {
    FCplot=FCplot+
      ggplot2::guides(color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))
  }


  png(filename=filename, width =1550, height = 1200, res=300)
  suppressWarnings(
    gridExtra::grid.arrange(FCplot, nrow=1, ncol=1,
                            left=grid::textGrob("Left hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = 0.5),
                            right=grid::textGrob("Right hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = -3))
  )
  dev.off()
}
########################################################################################################
########################################################################################################
#' @title edgelist
#' @description Listing the edges in a connectogram.
#'
#' @details This function takes a vector of edge-to-edge connections, and returns a data.frame object containing the edges (defined by their pair of connecting nodes) and their edge values
#'
#' @param data a vector of edge values with a length of 4005, 7021, 23871 or 30135.
#' @returns a data.frame() object with the columns of `node_1`, `node_2` and `weight` columns
#'
#' @examples
#'
#' results=sample(c(1,0, -1), 7021, replace = TRUE, prob = c(0.01, 0.98,0.01))
#' edgelist(data=results)
#'
#' @importFrom igraph graph_from_adjacency_matrix get.edgelist E
#' @export

########################################################################################################
########################################################################################################
edgelist=function(data)
{
  ##selecting atlas
  edge_lengths=c(4005,7021,23871,30135)


  if(is.na(match(length(data),edge_lengths)))
  {
    stop("The length of the input vector does not fit any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")
  } else
  {
    atlas=match(length(data),edge_lengths)
  }

  ##plot parameters
  label=labels_dat[[match(length(data),edge_lengths)]]
  label=label[order(label$oldorder),]
  label$labels=paste(label$hemi,"_",label$labels,sep="")
  nnodes=nrow(label)

  ##rehaping data into connectivity matrix

  conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
  conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = data
  conmat_NxN=conmat_NxNhalf+t(conmat_NxNhalf)

  nodeorder=as.numeric(rep(NA,nnodes))
  for (rowno in 1:nnodes){
    nodeorder[rowno]=which(label$neworder==rowno)
  }

  colnames(conmat_NxN)=label$labels
  rownames(conmat_NxN)=label$labels
  reordered=conmat_NxN[nodeorder,nodeorder]
  RegionsFC=as.factor(label$regionlabel)[nodeorder]

  ##graph object

  graphobj=igraph::graph_from_adjacency_matrix(reordered, mode="undirected", diag=FALSE, weighted=T)
  edgelist=data.frame(cbind(igraph::get.edgelist(graphobj) , igraph::E(graphobj)$weight))
  names(edgelist)=c("node_1","node_2","weight")
  return(edgelist)
}

########################################################################################################
########################################################################################################
