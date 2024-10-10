## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title vizChord
#'
#' @description Visualizing brain connectivity profiles with a chord diagram
#'
#' @details This function takes a matrix (NROW=number of edges in the connectome; NCOL=number of edges in the connectome) of edge values and visualizes the average network-to-network connectivity in a chord diagram.
#'
#' @param data a vector of edge values with a length of 78, 4005, 7021, 23871 or 30135.
#' @param width width (in pixels) of each connectogram. Set to 1200 by default. 
#' @param height height (in pixels) of each connectogram . Set to 1200 by default.
#' @param leg.height height (in pixels) of legend, in pixels. Set to 100 by default. Not used for single row data
#' @param ncol number of columns in the plot. Not used for single row data
#' @param nrow number of rows in the plot. Not used for single row data
#' @param title a vector of strings to be used as title
#' @param hot color for the positive connections.Set to `#F8766D` by default.
#' @param cold color for the negative connections.Set to `#00BFC4` by default.
#' @param colorscheme an optional vector of color names or color codes to color code the networks.
#' @param colorbar_title title for the colorbar legend
#' @param filename output filename with a *.png file extension. Set to `conn.png` by default
#'
#' @returns outputs a .png image
#'
#' @examples
#'
#' results=runif(7021, min = -1, max = 1)
#' vizChord(data=results, filename="FC_chord119.png")
#'
#' @importFrom ggplot2 ggplot aes scale_colour_gradient2 geom_point guides theme guide_colorbar element_text
#' @importFrom circlize colorRamp2 chordDiagram
#' @importFrom cowplot get_plot_component plot_grid
#' @importFrom ggplotify as.grob
#' @importFrom png readPNG
#' @importFrom grid grid.raster
#' @export
##Main function
########################################################################################################
########################################################################################################
vizChord=function(data, hot="#F8766D", cold="#00BFC4", width=1200, height=1200,filename="conn.png", colorscheme, title,leg.height=100, ncol=1,nrow=1, colorbar_title="Connectivity Strength")
{
  edge_lengths=c(4005,7021,23871,30135,78)
  if(NCOL(data)==1 | NROW(data)==1)
  {
    if(is.na(match(length(data),edge_lengths)))  
    {
      stop("The length of the input vector is not consisent with any of the recognized parcellation schemes. The input vector should contain 4005, 7021, 23871 or 30135 values")
    } else
    {
      atlas=match(length(data),edge_lengths)
      data=matrix(data,nrow=1)
    }  
  } else if(is.na(match(NCOL(data),edge_lengths)))
  {
    stop("The number of columns in the matrix is not consistent any of the recognized parcellation schemes. The input matrix should contain 4005, 7021, 23871 or 30135 columns")
  } else  
  {
    atlas=match(NCOL(data),edge_lengths)
    if((nrow*ncol)<NROW(data))
      {
       stop("Not enough columns or rows are specified to display the subplots")
      }
  } 
  
  if(missing("title")){title=rep(" ",NROW(data))}
  chordplots=list()
  
  if(NCOL(data)!=78)
  {
    if(missing("colorscheme") & NCOL(data)!=4005) {colorscheme = c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD")}
    else if(missing("colorscheme") & NCOL(data)==4005)  {colorscheme = c("#D53E4F","#FC8D59","#FEE08B","#FFFFBF","#E6F598","#99D594","#3288BD")}
    
    for(rowno in 1:NROW(data))  {chordplots[[rowno]]=genChord(data[rowno,],hot,cold,colorscheme)}
    
  } else
  {
    #set default colors
    if(missing("colorscheme"))  {colorscheme=c("#A71B4BFF","#D04939FF","#EB7803FF","#F5A736FF","#FBCF6FFF","#FEF1A6FF","#E2F8B5FF","#9CE5ADFF","#43CBB1FF","#00AAB6FF","#0080B2FF","#584B9FFF")}
    
    for(rowno in 1:NROW(data))  {chordplots[[rowno]]=genChord_12x12(data[rowno,],hot,cold,colorscheme)}
  }
  ##make colorbar legend
  leg.dat=data.frame(v1=-100:100)
  legend.plot=ggplot2::ggplot(leg.dat, ggplot2::aes(color=v1, x=v1, y=v1))+
              ggplot2::scale_colour_gradient2(name=colorbar_title,low=cold,mid="white",high=hot,
                                      guide = "colourbar", limits=c(-1,1), breaks=c(-1,1),
                                      labels=c("Strong negative","Strong positive"))+
              ggplot2::geom_point()+
              ggplot2::guides(color=ggplot2::guide_colorbar(ticks=F,title.position = "left",barheight = 0.5))+
              ggplot2::theme(legend.position = "bottom",legend.title=ggplot2::element_text(face="bold", size=7,vjust=1),
                             legend.text =ggplot2::element_text(size=6))
  legend = cowplot::get_plot_component(legend.plot, 'guide-box-bottom', return_all = TRUE)
    
  ##combine chord plots and legend  
  chord=do.call(eval(parse(text="cowplot::plot_grid")), args = list(plotlist=chordplots,nrow=nrow, ncol=ncol,labels=title))
  main=cowplot::plot_grid(chord,legend,ncol=1, nrow=2, rel_heights = c(nrow*height,leg.height))
  
  png(filename=filename, width=ncol*width,height=nrow*height+leg.height,res=300)
    print(main)
  dev.off()
  img=png::readPNG(source =filename)
  grid::grid.raster(img)
}
########################################################################################################
########################################################################################################
genChord_12x12=function(data,hot,cold, colorscheme)
{
  localenv = environment()
  
  nodesname<-c("Auditory","Cingulo-\nopercular","Cingulo-\nparietal",
               "Default\nmode","Dorsal\nattention","Frontoparietal",
               "Retrosplenial\ntemporal","Sensorimotor\nhand","Sensorimotor\nmouth",
               "Salience","Ventral\nattention","Visual")
  no_nodes=12
  
  #color parameters
  colarrFC=array()
  pos_color_range= colorRampPalette(c("white",hot))
  neg_color_range = colorRampPalette(c("white",cold))
  pos_color_val=pos_color_range(101)
  neg_color_val=neg_color_range(101)
  
  #reconstruct 12 x 12 FC matrices
  FC_12X12=array(rep(NA,12^2),dim=c(12,12))
  FC_12X12[upper.tri(FC_12X12, diag=T)] = data
  
  a=1
  
  datFC0<-array(dim=c(144,3))
  for (rowno in 1:no_nodes)
  {
    for (colno in 1:no_nodes)
    {
      datFC0[a,1]=rowno
      datFC0[a,2]=colno
      datFC0[a,3]=FC_12X12[rowno,colno]
      a=a+1
    }
  }
  #converting 12 x 12 FC matrices to a 'from node to node' data frame
  datFC=datFC0[-which(is.na(datFC0[,3])),]
  datFC=data.frame(datFC, stringsAsFactors = F)
  
  colnames(datFC)<-c("from","to","value")
  
  #sets up the color function such that darker colors correspond to stronger connections
  
  col_funFC = circlize::colorRamp2(range(datFC$value), c(hot,cold))
  
  for (label in 1:no_nodes){
    datFC$from[which(datFC$from==label)]=nodesname[label]
    datFC$to[which(datFC$to==label)]=nodesname[label]
  }
  
  for (datval in which(datFC$value>0))  {colarrFC[datval]<-pos_color_range[round(abs(datFC$value[datval])/max(abs(datFC$value))*100)+1]}
  for (datval in which(datFC$value<0))  {colarrFC[datval]<-neg_color_range[round(abs(datFC$value[datval])/max(abs(datFC$value))*100)+1]}
  
  FCchord=ggplotify::as.grob(~circlize::chordDiagram(datFC, col=colarrFC, self.link = 1, grid.col=colorscheme[1:12], annotationTrack = c("grid","name"),link.border=colarrFC),envir = localenv)
  return(FCchord)
}

########################################################################################################
########################################################################################################

genChord=function(data,hot,cold, colorscheme)
{
  localenv <- new.env()
  
  edge_lengths=c(4005,7021,23871,30135)
  labels_dat=get('labels_dat')
  label=labels_dat[[match(length(data),edge_lengths)]]
  label=label[order(label$region),]
  labelnameFC=unique(label$regionlabel)
  
  regionnoFC=label$region
  noregions=length(labelnameFC)
  nnode=nrow(label)
  
  ##color parameters
  pos_color_range= colorRampPalette(c("white",hot))
  neg_color_range = colorRampPalette(c("white",cold))
  pos_color_val=pos_color_range(101)
  neg_color_val=neg_color_range(101)
  
  #if there are more colors in the colorscheme that the number of regions
  if(noregions>length(colorscheme)) {colorscheme=colorscheme[1:noregions]}
  
  FC_matrix=array(rep(0,nnode^2),dim=c(nnode,nnode))
  FC_matrix[upper.tri(FC_matrix, diag=FALSE)] = data
  FC_matrix=FC_matrix+t(FC_matrix)
  
  FCregionmat=array(dim=c(noregions,noregions))
  for (rowno in 1:noregions)
  {
    for (colno in 1:noregions)  {FCregionmat[rowno,colno]=mean(FC_matrix[label$oldorder[which(regionnoFC==rowno)],label$oldorder[which(regionnoFC==colno)]])}
  }
  
  #region x region matrix FC
  FCregionmat[lower.tri(FCregionmat)]=NA
  x=array(dim=c(noregions^2,3))
  
  a=1
  for (rowno in 1:noregions)
  {
    for (colno in 1:noregions)
    {
      x[a,1]=rowno
      x[a,2]=colno
      if (rowno==colno) {x[a,3]=(FCregionmat[rowno,colno])/2}
      else {x[a,3]=FCregionmat[rowno,colno]}
      a=a+1
    }
  }
  
  datFC=x[-which(is.na(x[,3])),]
  datFC=data.frame(datFC, stringsAsFactors = F)
  
  colnames(datFC)=c("from","to","value")
  col_funFC = circlize::colorRamp2(range(datFC$value), c(cold,hot))
  
  for (label in 1:noregions){
    datFC$from[which(datFC$from==label)]=labelnameFC[label]
    datFC$to[which(datFC$to==label)]=labelnameFC[label]
  }
  colarrFC=array()
  
  if(NROW(which(datFC$value==0))==0){datrmFC=datFC} 
  else {datrmFC=datFC[-which(datFC$value==0),]}
  
  if(length(unique(c(datrmFC$from,datrmFC$to)))< noregions)  {stop(paste("Connections are absent in",noregions-length(unique(c(datrmFC$from,datrmFC$to))),"of the networks. Please use a connectogram instead",sep=" "))}
  
  for (datval in which(datrmFC$value>0))  {colarrFC[datval]=pos_color_val[round(abs(datrmFC$value[datval])/max(abs(datrmFC$value))*100)+1]}
  for (datval in which(datrmFC$value<0))  {colarrFC[datval]=neg_color_val[round(abs(datrmFC$value[datval])/max(abs(datrmFC$value))*100)+1]}
  
  FCchord=ggplotify::as.grob(~suppressWarnings(circlize::chordDiagram(datrmFC, col=colarrFC, order=labelnameFC, self.link = 2,grid.col=colorscheme,annotationTrack = c("grid","name"),link.border=colarrFC,)),envir = localenv)

  return(FCchord)
}
