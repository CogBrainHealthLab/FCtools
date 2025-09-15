## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
#' @title vizGlassbrain
#'
#' @description Visualizing brain connectivity with a glass brain plot
#'
#' @details This function takes a vector of edge values and visualizes the edge-to-edge connectivity with a glass brain plot
#'
#' @param data a vector of edge values with a length of 4005, 7021, 23871 or 30135
#' @param surf_color  color of the cortical surface. Set to `'grey'` by default
#' @param node_color  color of the nodes. Set to `'#00BA38'` by default
#' @param surf_alpha  alpha value of the cortical surface, where 0 will cause the cortical surface to disappear and 1 will cause the cortical surface to be completely opaque. Set to `0.2` by default
#' @param cmap A string vector containing 2 to 3 color names/codes specifying the colors to be used for the color scale. See `RColorBrewer::display.brewer.all()` for all possible cmap options. If none are specified, appropriate colors will be automatically selected according to `range(data)`
#' @param node_size size parameter for the dots representing the nodes. Set to `8` by default.
#' @param node_label option to show node labels. Set to `TRUE` by default.
#' @param node_label_size font size of the node labels. Set to `10` by default.
#' @param node_label_color font color of the node labels. Set to `black` by default.
#' @param edgethickness a value to adjust the thickness of the edges. Set to `8` by default.
#' @param limits A combined pair of numeric values composed of the lower and upper color scale limits of the plot. When left unspecified, the symmetrical limits `c(-max(abs(data),max(abs(data)))` will be used. 
#' @param colorbar_title title for the colorbar legend. Set to `'Connectivity strength` by default
#' @param orientation_labels A boolean object specifying if orientation labels are to be displayed. Set to `TRUE` by default
#' @param remove_brain A boolean object specifying cortical surface should be removed. Set to `FALSE` by default
#'
#' @returns outputs a plot_ly object
#'
#' @examples
#' 
#' mask=sample(c(1,0), 7021, replace = TRUE, prob = c(0.001, 0.999))
#' data=runif(7021,min = -1,max=1)*mask
#'
#' vizGlassbrain(data,orientation_labels = T)
#' 
#' @importFrom plotly plot_ly add_trace layout
#' @export
########################################################################################################
########################################################################################################
vizGlassbrain=function(data,
                       surf_color="grey",
                       node_color="#00BA38",
                       surf_alpha=0.2,
                       cmap,
                       node_size=8,
                       node_label=TRUE,
                       node_label_size=10,
                       node_label_color="black",
                       edgethickness=8,
                       limits, 
                       colorbar_title="Connectivity strength",
                       orientation_labels=TRUE,
                       remove_brain=FALSE)
{
##color scale parameters
  #defining cmap if missing
  if(missing("cmap"))
  {
    if(range(data,na.rm = TRUE)[1]>=0)  {cmap=c("#F5FACD","#F8766D")}
    else if (range(data,na.rm = TRUE)[2]<=0)  {cmap=c("#619CFF","#E7F1D5")}
    else  {cmap=c("#619CFF","white","#F8766D")}
  }
  #enabling custom color scales
  if(length(cmap)==2) {cmap=list(list(0,cmap[1]), list(1,cmap[2]))
  } else if (length(cmap)==3) {cmap=list(list(0,cmap[1]), list(0.5,cmap[2]),list(1,cmap[3]))}
  
  #setting color scale limits
  maxlimit=max(abs(range(data,na.rm = TRUE)))
  if(missing(limits)) 
  {
    limits.range=range(data,na.rm = T)
    if(limits.range[1]>=0) {limits=c(0,limits.range[2])} ##if image contains all positive values
    else if(limits.range[2]<=0) {limits=c(limits.range[1],0)} ##if image contains all negative values
    else if(limits.range[1]<0 & limits.range[2]>0){limits=c(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values
  } else {limits=c(limits[1],limits[2])}

##reshaping FC vector to FC matrix
  labels_dat=get('labels_dat')
  lengths=c(4005,7021,23871,30135)
  nodes=c(90,119,219,246)
  if(is.na(match(length(data),lengths)))  {stop("The length of the input vector is not consisent with any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")}
  
  atlasno=match(length(data),lengths)
  atlas=labels_dat[[atlasno]]
  
  FCmat=matrix(NA,nodes[atlasno],nodes[atlasno])
  FCmat[upper.tri(FCmat)]=data
  FCmat[lower.tri(FCmat)] <- t(FCmat)[lower.tri(FCmat)]
  
  FCmat.bin=FCmat
  FCmat.bin[FCmat.bin!=0]=1

##extracting MNI coords  
  edgelist=igraph::as_edgelist(igraph::graph_from_adjacency_matrix(FCmat.bin))
  MNIdat=matrix(NA,nrow=NROW(edgelist)*2,ncol=5)
  labs=rep(NA,NROW(edgelist)*2)
  count=1
  for(edge in 1:NROW(edgelist))
  {
    MNIdat[count,1]=edge
    MNIdat[count,2:4]=as.numeric(atlas[edgelist[edge,1],7:9])
    MNIdat[count,5]=FCmat[edgelist[edge,][1],edgelist[edge,][2]]
    labs[count]=atlas$labels[edgelist[edge,1]]
    count=count+1
    
    MNIdat[count,1]=edge
    MNIdat[count,2:4]=as.numeric(atlas[edgelist[edge,2],7:9])
    MNIdat[count,5]=FCmat[edgelist[edge,][1],edgelist[edge,][2]]
    labs[count]=atlas$labels[edgelist[edge,2]]
    count=count+1
  }
  MNIdat=data.frame(MNIdat)
  colnames(MNIdat)=c("edge","X","Y","Z","strength")  
  
  ##create blank cortical surface
  fs5brain=get('fs5brain')
  fig=plotly::plot_ly()
  
  if(remove_brain==F)
  {
  fig=plotly::add_trace(fig,type = 'mesh3d',hoverinfo="skip",
                        x = fs5brain[[1]][,1],
                        y = fs5brain[[1]][,2],
                        z = fs5brain[[1]][,3],
                        i = fs5brain[[2]][, 1] - 1,  # plotly uses 0-based indexing, so subtract 1
                        j = fs5brain[[2]][, 2] - 1,
                        k = fs5brain[[2]][, 3] - 1,facecolor=rep(surf_color,NROW(fs5brain[[2]])),opacity=surf_alpha)
  }
  ##add nodes and edges  
  fig=plotly::add_trace(fig,data = MNIdat, x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines+markers', split = ~edge,
                        line = list(width = edgethickness, 
                               color = ~strength, 
                               colorscale = cmap,
                               cauto = F,
                               cmin = limits[1],
                               cmax = limits[2],
                               colorbar = list(ticklabelposition="outside left",
                                               title = list(text=colorbar_title,side="right"))),
                        text=labs,hoverinfo = 'text',
                        marker = list(size = node_size, color = node_color),showlegend = F)

  ##add node labels
  if(node_label==T)
  {
  fig=plotly::add_text(fig,data = MNIdat, x = ~X, y = ~Y, z = ~Z, text = labs, textfont=list(size=node_label_size, color=node_label_color))
  }
  ##turn off grid, tick labels, and axes labels
  fig=plotly::layout(fig,hoverlabel = list(align = "left"),
                     scene = list(camera=list(eye = list(x = 0, y = 1.5, z = 1.5)),
                                  xaxis = list(showgrid = F,showticklabels=F,showspikes=F,zeroline=F, title=""),
                                  yaxis = list(showgrid = F,showticklabels=F,showspikes=F,zeroline=F, title=""),
                                  zaxis = list(showgrid = F,showticklabels=F,showspikes=F,zeroline=F, title="")))
  
  ##add optional orientation labels
  if(orientation_labels==T)
  {
    fig=plotly::layout(fig,hoverlabel = list(align = "left"),
                       scene = list(xaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title=""),
                                    yaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title=""),
                                    zaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title="")))
    axx = list(ticketmode = 'array',ticktext =  c("Left","Right"),tickvals = range(fs5brain[[1]][,1]))
    axy = list(ticketmode = 'array',ticktext = c("Posterior","Anterior"),tickvals = range(fs5brain[[1]][,2]))
    axz = list(ticketmode = 'array',ticktext = c("Inferior","Superior"),tickvals = range(fs5brain[[1]][,3]))
    
    fig = plotly::layout(fig,scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  }
  return(fig)
}
