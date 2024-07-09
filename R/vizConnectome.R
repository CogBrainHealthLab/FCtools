## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
#' @title vizConnectogram
#'
#' @description Visualizing mbrain connectivity profiles with multiple connectogram
#'
#' @details This function takes a matrix (NROW=number of connectogram; NCOL=number of edges in the connectogram) of edge values and visualizes the edge-to-edge connectivity with multiple connectograms
#'
#' @param data a matrix of edge values with 4005, 7021, 23871 or 30135 columns
#' @param ncol number of columns in the plot. Not used for single row data
#' @param nrow number of rows in the plot. Not used for single row data
#' @param hot color for the positive connections. Set to `#F8766D` by default.
#' @param cold color for the negative connections. Set to `#00BFC4` by default.
#' @param colorscheme an optional vector of color names or color codes to color code the networks.
#' @param filename output filename with a *.png file extension. Set to "conn.png" by default
#' @param edgethickness a value to adjust the thickness of the edges. Set to 0.8 by default.
#' @param title a vector of strings to be used as title
#' @param width width (in pixels) of each connectogram. Set to 1000 by default. Not used for single row data
#' @param height height (in pixels) of each connectogram . Set to 1050 by default. Not used for single row data
#' @param leg.height height (in pixels) of legend, in pixels. Set to 100 by default. Not used for single row data
#' @param margins vector of 4 values specifying the amount of empty space on the left, right, top and bottom for each connectogram. Set to `c(1.2,1.2,1,1.5)` by default. You might want to adjust these values if the text labels get cut off by a neighbouring connectogram or the legend. Not used for single row data
#' @returns outputs a .png image
#'
#' @examples
#' if (FALSE) {
#' data=matrix(sample(c(1,0, -1), 23871*6, replace = TRUE, prob = c(0.001, 0.998,0.001)),nrow=6)
#' vizConnectogram(data=data,ncol=3, nrow=2)
#' }
#' @importFrom ggplot2 ggplot aes scale_colour_gradient2 geom_point guides theme guide_colorbar element_text scale_color_manual expand_limits unit guide_legend
#' @importFrom igraph graph_from_adjacency_matrix get.edgelist E edge_attr
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar
#' @importFrom ggraph ggraph geom_edge_arc scale_edge_alpha_continuous scale_edge_color_manual geom_node_point geom_node_text theme_graph
#' @importFrom cowplot plot_grid get_plot_component
#' @export

########################################################################################################
########################################################################################################
  vizConnectogram=function(data,
                           hot="#F8766D", 
                           cold="#00BFC4",
                           ncol,
                           nrow, 
                           edgethickness=0.8,
                           filename="conn.png", 
                           colorscheme, 
                           title,
                           width=1000,
                           height=1050,
                           leg.height=200,
                           margins=c(1.2,1.2,1,1.5))
{
  ##selecting atlas
  edge_lengths=c(4005,7021,23871,30135)
  if(NCOL(data)==1 | NROW(data)==1)
  {
    if(is.na(match(length(data),edge_lengths)))
    {
      stop("The length of the input vector is not consisent with any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")
    } else
    {
    atlas=match(length(data),edge_lengths)
    mode="vector"
    data=matrix(data,nrow=1)
    }  
  } else if(is.na(match(NCOL(data),edge_lengths)))
  {
    stop("The number of columns in the matrix is not consistent any of the recognized parcellation schemes. The input matrix should contain 4005, 7021, 23871 or 30135 columns")
  } else
  {
    atlas=match(NCOL(data),edge_lengths)
    mode="matrix"
  } 
  
  ##parameters for all atlases
  param=list()
  param$nodecol=list(c("#D53E4F","#FC8D59","#FEE08B","#FFFFBF","#E6F598","#99D594","#3288BD"),
                     c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD"))
  param$nodecol[3:4]=param$nodecol[2]
  param$nodesize=c(1.4,1.3,1,1)
  param$margins=c(-margins[1], margins[2],margins[3],-margins[4])
  label=labels_dat[[atlas]]
  param$nodelevels=unique(label$regionlabel)
  label=label[order(label$oldorder),]
  label.neworder=label[order(label$neworder),]
  param$nodelevels=unique(label.neworder$regionlabel)
  label$regionlabel = factor(label$regionlabel,levels = param$nodelevels)
  if(missing("colorscheme")){colorscheme = param$nodecol[[atlas]]}
  if(missing("title")){title=rep(" ",NROW(data))}
  param$xlim=list(c(-1.2, 1.2),
                  c(-1.1, 1.1),
                  c(-1.1, 1.1),
                  c(-1, 0.7))
  param$ylim=list(c(-1.2, 1.3),
                  c(-1.2, 1.1),
                  c(-1.2, 1.1),
                  c(-1, 1.1))
  nnodes=NROW(label)
  
  ##generate first plot
    #reshaping FC vector to FC matrix
    conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
    conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = data[1,]
    conmat_NxN=conmat_NxNhalf+t(conmat_NxNhalf)
    
    nodeorder=as.numeric(rep(NA,nnodes))
    for (rowno in 1:nnodes) {nodeorder[rowno]=which(label$neworder==rowno)}
    
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
    
  ##different steps depending on mode
  if(mode=="matrix")
  {
    legend.plot=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
      ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness*1.5, show.legend = T,) +
      ggraph::scale_edge_alpha_continuous(guide="none")+
      ggraph::scale_edge_color_manual(name="Edges", labels=c("Positive","Negative"),values=c(hot,cold))+
      ggplot2::scale_color_manual(values =colorscheme, name="Network")+
      ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=param$nodesize[atlas]*1.5, shape=19,show.legend = T) +
      ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                          angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                          hjust = ifelse(x > 0, 0 ,1)), size=param$nodesize[atlas]) +
      ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
      ggplot2::guides(edge_color = ggplot2::guide_legend(override.aes = list(shape = NA),nrow=2),
                      color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))+
      ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4),
                     legend.title=ggplot2::element_text(size=8,face = "bold",hjust=0.5),
                     legend.text=ggplot2::element_text(size=7),
                     legend.position ="bottom", legend.title.position = "top",
                     legend.key.height = ggplot2::unit(c(0, 0, 0, 0), "cm"))
    
    legend = cowplot::get_plot_component(legend.plot, 'guide-box-bottom', return_all = TRUE)
    
    ##generating actual plots
    ggplot.obj=list()
    for(row_dat in 1:NROW(data))
    {
      ggplot.obj[[row_dat]]=suppressWarnings(genplot(row_data = data[row_dat,],
                                                     title=title[row_dat],
                                                     nnodes=nnodes,
                                                     label=label,
                                                     edgethickness=edgethickness, 
                                                     hot=hot, 
                                                     cold=cold,
                                                     colorscheme=colorscheme,
                                                     param=param,
                                                     atlas=atlas))
    }
    
    main=gridExtra::grid.arrange(grobs=ggplot.obj,nrow=nrow,ncol=ncol)
    
    png(filename=filename,width=ncol*width, height=(nrow*height)+leg.height,res=300)
    suppressWarnings(print(cowplot::plot_grid(main,legend,nrow=2, rel_heights = c((nrow*height),leg.height))))
    dev.off()
  } else if(mode=="vector")
  {
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
                     legend.position.inside = c(1.1,0),legend.justification=c(1,0), legend.key.height = ggplot2::unit(c(0, 0, 0, 0), "cm"))+
      ggplot2::guides(edge_color = ggplot2::guide_legend(override.aes = list(shape = NA)),color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))
    
      png(filename=filename, width =1600, height = 1200, res=300)
      suppressWarnings(
        gridExtra::grid.arrange(FCplot, nrow=1, ncol=1,
                                left=grid::textGrob("Left hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = 0.5),
                                right=grid::textGrob("Right hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = -5))
      )
    dev.off()
  }
}
##function to generate individual plots
genplot=function(row_data,title,nnodes,label,edgethickness, hot, cold,colorscheme,param, atlas)
{
  conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
  conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = row_data
  conmat_NxN=conmat_NxNhalf+t(conmat_NxNhalf)
  
  nodeorder=as.numeric(rep(NA,nnodes))
  for (rowno in 1:nnodes) {nodeorder[rowno]=which(label$neworder==rowno)}
  
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
  ggplot.obj=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
    ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness, show.legend = F) +
    ggraph::scale_edge_alpha_continuous(guide="none")+
    ggraph::scale_edge_color_manual(name="Edges", labels=c("Positive","Negative"),values=c(hot,cold))+
    ggplot2::scale_color_manual(values =colorscheme, name="Network")+
    ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=param$nodesize[atlas], shape=19,show.legend = F) +
    ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                        angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                        hjust = ifelse(x > 0, 0 ,1)), size=param$nodesize[atlas]) +
    ggplot2::ggtitle(title)+
    ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
    ggplot2::expand_limits(x = param$margin[1:2], y = param$margin[3:4])+
    ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4))
}
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
