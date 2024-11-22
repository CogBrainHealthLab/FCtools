## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
#' @title vizConnectogram
#'
#' @description Visualizing brain connectivity profiles with multiple connectogram
#'
#' @details This function takes a matrix (NROW=number of edges in the connectome; NCOL=number of edges in the connectome) of edge values and visualizes the edge-to-edge connectivity with multiple connectograms
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
#' @param node.size size parameter for the dots representing the nodes. If not specified, an appropriate size will be set.
#' @param node.text.size size parameter for the text labels on the nodes. Set to 1 by default.
#' @param title.size size parameter for the title. Set to 11 by default.
#' @param title.alignment string object specifying title's alignment relatively to its plot options are 'left', 'center' (default), 'right'.
#' @param legend.title.size size parameter for the legend title. Set to 8 by default.
#' @param legend.text.size size parameter for the legend text. Set to 6 by default.
#' @param limits a pair of values that governs the limits of the edge strengths displayed. If missing `limits=range(abs(data),na.rm=T)`
#' @param colorbar_title title for the colorbar legend
#' @param edge_labels Vector of two strings defining the labels for the edge legends. Default is c("Positive","Negative").
#' @param row_title a vector of strings to be used as left vertical titles for each row of plots when there are many
#' @returns outputs a .png image
#'
#' @examples
#' if (FALSE) {
#' data=matrix(sample(c(1,0, -1), 23871*6, replace = TRUE, prob = c(0.001, 0.998,0.001)),nrow=6)
#' vizConnectogram(data=data,ncol=3, nrow=2)
#' }
#' @importFrom ggplot2 ggplot aes scale_colour_gradient2 geom_point guides theme guide_colorbar element_text scale_color_manual expand_limits unit guide_legend coord_fixed
#' @importFrom igraph graph_from_adjacency_matrix get.edgelist E edge_attr
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar grid.raster nullGrob textGrob
#' @importFrom ggraph ggraph geom_edge_arc scale_edge_alpha_continuous scale_edge_color_manual geom_node_point geom_node_text theme_graph
#' @importFrom cowplot plot_grid get_plot_component
#' @importFrom png readPNG
#' @export

########################################################################################################
########################################################################################################
vizConnectogram=function(data,
                         hot="#F8766D", 
                         cold="#00BFC4",
                         ncol=1,
                         nrow=1, 
                         edgethickness=0.8,
                         filename="conn.png", 
                         colorscheme, 
                         title,
                         width=1000,
                         height=1050,
                         leg.height=200,
                         margins=c(1.2,1.2,1,1.5),
                         node.size,
                         node.text.size=1,
                         legend.text.size=7,
                         legend.title.size=8,
                         limits,
                         title.size=11,
                         title.alignment='center',
                         colorbar_title="Edge Strength",
                         edge_labels=c("Positive","Negative"),
                         row_title)
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
  data("labels_dat", package = "FCtools")
  label=labels_dat[[atlas]]
  label=label[order(label$oldorder),]
  label.neworder=label[order(label$neworder),]
  param$nodelevels=unique(label.neworder$regionlabel)
  label$regionlabel = factor(label$regionlabel,levels = param$nodelevels)
  if(missing("colorscheme")){colorscheme = param$nodecol[[atlas]]}
  if(missing("title")){title=rep(" ",NROW(data))}
  if(missing("node.size")){node.size=param$nodesize[atlas]}
  if(missing("limits")){limits=range(abs(data),na.rm=T)}
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
  
  edgelab=edge_labels
  if(length(which(posnegFC=="1_pos"))==0)
  {
    hot=cold
    edgelab=edge_labels[2]
  }
  
  igraph::edge_attr(graphobjFC, "weight", index = igraph::E(graphobjFC))=abs(EvalFC)
  igraph::edge_attr(graphobjFC, "posFC", index = igraph::E(graphobjFC))=posnegFC
  
  ##different steps depending on mode
  if(mode=="matrix")
  {
    if((nrow*ncol)<NROW(data))  {stop("Not enough columns or rows are specified to display the subplots")}
    if(range(data)[1]<0 & range(data)[2]>0)
    {
      #generate dummy dat for legend
      #reshaping FC vector to FC matrix
      edgelab=edge_labels
      conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
      conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = sample(c(-1,1),size=length(data[1,]),replace=T,prob = c(50,50))
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
    }
    
    legend.plot=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
      ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness*1.5, show.legend = TRUE) +
      ggraph::scale_edge_alpha_continuous(guide="none", limits=limits)+
      ggraph::scale_edge_color_manual(name=colorbar_title, labels= edgelab,values=c(hot,cold),drop = FALSE)+
      ggplot2::scale_color_manual(values =colorscheme, name="Network")+
      ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=node.size*1.5, shape=19,show.legend = T) +
      ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                          angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                          hjust = ifelse(x > 0, 0 ,1)), size=node.text.size) +
      ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
      ggplot2::guides(edge_color = ggplot2::guide_legend(override.aes = list(shape = NA),nrow=2),
                      color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))+
      ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4), 
                     legend.title=ggplot2::element_text(size=legend.title.size,face = "bold",hjust=0.5),
                     legend.text=ggplot2::element_text(size=legend.text.size),
                     legend.position ="bottom", legend.title.position = "top",
                     legend.key.height = ggplot2::unit(c(0, 0, 0, 0), "cm"))
    
    legend = suppressWarnings(cowplot::get_plot_component(legend.plot, 'guide-box-bottom', return_all = TRUE))
    
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
                                                     atlas=atlas,
                                                     node.text.size=node.text.size,
                                                     node.size=node.size,
                                                     title.size=title.size,
                                                     title.alignment=title.alignment,
                                                     colorbar_title=colorbar_title,
                                                     edge_labels=edge_labels,
                                                     limits=limits
                                                     ))
                                                     
      }
    
    #main=suppressWarnings()
    
    #If row titles specified:
    if (!missing(row_title)) 
    {
      row_titles=list()
      #Create rotated side title for each row title provided
      for (rowt in row_title) 
      {row_titles <- append(row_titles, 
          list(grid::textGrob(rowt, rot = 90,
                              gp = grid::gpar(fontsize = title.size))))
      }

      #keep only as many titles as there are rows
      grobs <- list()
      for (i in 1:nrow) 
      {if (i <= length(row_titles))
        {grobs <- c(grobs, list(row_titles[[i]]))} 
      else 
        {grobs <- c(grobs, list(grid::textGrob("")))} #if not as many
      
          #depending on how many plots per row (user-chosen ncol), 
          #add title and then plots for each row
          for (j in 1:ncol) { 
            index <- (i - 1) * ncol + j;
            if (index <= length(ggplot.obj)) 
            {grobs <- c(grobs, list(ggplot.obj[[index]]))} 
            #if missing plot
            else {grobs <- c(grobs, list(grid::nullGrob()))} 
            }
      }
      ncol=ncol+1 #first column will be dedicated to title
      width=(ncol-1)*width
      widths = c(0.2, rep(5, ncol-1))
      ggplot.obj=grobs
    }
    else { 
      #Default width if no row title
      width=ncol*width
      widths=rep(1,ncol)}
    
    png(filename=filename,width=width, height=(nrow*height)+leg.height,res=300)
    suppressWarnings(print(cowplot::plot_grid(gridExtra::grid.arrange(grobs=ggplot.obj,nrow=nrow,ncol=ncol, widths = widths),legend,nrow=2, rel_heights = c((nrow*height),leg.height))))
    dev.off()
    img=png::readPNG(source =filename)
    grid::grid.raster(img)
  } else if(mode=="vector")
  {
    FCplot=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
      ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness, show.legend = T) +
      ggraph::scale_edge_alpha_continuous(guide="none", limits=limits)+
      ggraph::scale_edge_color_manual(name=colorbar_title, labels=edgelab,values=c(hot,cold))+
      ggplot2::scale_color_manual(values =colorscheme, name="Network")+
      ggplot2::ggtitle(title[1])+
      ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=node.size, shape=19,show.legend = T) +
      ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                          angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                          hjust = ifelse(x > 0, 0 ,1)), size=node.text.size) +
      ggplot2::coord_fixed()+
      ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
      ggplot2::expand_limits(x = param$xlim[[atlas]], y = param$ylim[[atlas]])+
      ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4),plot.title=ggplot2::element_text(size=title.size),
                     legend.title=ggplot2::element_text(size=legend.title.size-2,face = "bold"),
                     legend.text=ggplot2::element_text(size=legend.text.size-2),
                     legend.position.inside = c(1.1,0),legend.justification=c(1,0), legend.key.height = ggplot2::unit(c(0, 0, 0, 0), "cm"))+
      ggplot2::guides(edge_color = ggplot2::guide_legend(override.aes = list(shape = NA)),color= ggplot2::guide_legend(override.aes = list(edge_width = NA)))
    
    png(filename=filename, width =1600, height = 1200, res=300)
    suppressWarnings(
      gridExtra::grid.arrange(FCplot, nrow=1, ncol=1,
                              left=grid::textGrob("Left hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = 0.5),
                              right=grid::textGrob("Right hemisphere",gp = grid::gpar(fontface=2,fontsize = 6),rot=90, hjust = 0.5,x = -5))
    )
    dev.off()
    img=png::readPNG(source =filename)
    grid::grid.raster(img)
  }
}
##function to generate individual plots
genplot=function(row_data,
                 title,
                 nnodes,
                 label,
                 edgethickness, 
                 hot, 
                 cold,
                 colorscheme,
                 param, 
                 atlas,
                 node.text.size,
                 node.size,
                 title.size,
                 title.alignment,
                 colorbar_title,
                 edge_labels,
                 limits)
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
  
  if(length(which(posnegFC=="1_pos"))==0)
  {
    hot=cold
  }
   
  igraph::edge_attr(graphobjFC, "weight", index = igraph::E(graphobjFC))=abs(EvalFC)
  igraph::edge_attr(graphobjFC, "posFC", index = igraph::E(graphobjFC))=posnegFC
  
  #plot
  
  #alignment value
  if (title.alignment=='center') {talign=0.5}
  else if (title.alignment=='left') {talign=0}
  else if (title.alignment=='right') {talign=1}
  else {talign=0.5}
  
  ggplot.obj=ggraph::ggraph(graphobjFC, layout = 'linear', circular = TRUE) +
    ggraph::geom_edge_arc(ggplot2::aes(color=posnegFC, alpha=weight), edge_width=edgethickness, show.legend = F) +
    ggraph::scale_edge_alpha_continuous(guide="none",limits=limits)+
    ggraph::scale_edge_color_manual(name=colorbar_title, labels=edge_labels,values=c(hot,cold))+
    ggplot2::scale_color_manual(values =colorscheme, name="Network")+
    ggraph::geom_node_point(ggplot2::aes(colour = RegionsFC),size=node.size, shape=19,show.legend = F) +
    ggraph::geom_node_text(ggplot2::aes(label = name, x = x * 1.03, y = y* 1.03,
                                        angle = ifelse(atan(-(x/y))*(180/pi) < 0,90 + atan(-(x/y))*(180/pi), 270 + atan(-x/y)*(180/pi)),
                                        hjust = ifelse(x > 0, 0 ,1)), size=node.text.size) +
    ggplot2::coord_fixed()+
    ggplot2::ggtitle(title)+
    ggraph::theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
    ggplot2::expand_limits(x = param$margin[1:2], y = param$margin[3:4])+
    ggplot2::theme(plot.margin = rep(ggplot2::unit(0,"null"),4),
                   plot.title=ggplot2::element_text(size=title.size,
                                                    hjust=talign))
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
