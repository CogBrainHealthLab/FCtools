## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title vizChord
#'
#' @description Visualizing brain connectivity profiles with a chord diagram
#'
#' @details This function takes a vector of edges and visualizes the average network-to-network connectivity in a chord diagram.
#'
#' @param data a vector of edge values with a length of 78, 4005, 7021, 23871 or 30135.
#' @param width width of the chord diagram.Set to 1800 by default.
#' @param height height of the chord diagram.Set to 1800 by default.
#' @param hot color for the positive connections.Set to `#F8766D` by default.
#' @param cold color for the negative connections.Set to `#00BFC4` by default.
#' @param colorscheme an optional vector of color names or color codes to color code the networks.
#' @param filename output filename with a *.png file extension. Set to `conn.png` by default
#' @param colorbar `TRUE` or `FALSE` option on whether a colorbar for the connectivity strength will be displayed. Set to `TRUE` by default.
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
#' @importFrom cowplot get_legend plot_grid
#' @importFrom ggplotify as.grob
#' @export

##Main function
########################################################################################################
########################################################################################################
vizChord=function(data, hot="#F8766D", cold="#00BFC4", width=1800, height=1800,filename="conn.png", colorscheme, colorbar=TRUE)
{
  ##redirect to other functions depending on length of data
  if((length(data)==23871 |length(data)==7021)| length(data)==30135)
  {
    if(missing("colorscheme")){colorscheme = c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD")}

    vizChord_edge(data=data,hot=hot,cold=cold,width=width,height=height,filename=filename, colorscheme = colorscheme, colorbar=colorbar)

  } else if (length(data)==4005)
  {
    if(missing("colorscheme")){colorscheme = c("#D53E4F","#FC8D59","#FEE08B","#FFFFBF","#E6F598","#99D594","#3288BD")}

    vizChord_edge(data=data,hot=hot,cold=cold,width=width,height=height,filename=filename, colorscheme = colorscheme, colorbar=colorbar)
  } else if(length(data)==78)
  {
    if(missing("colorscheme"))
    {
      colorscheme=c("#A71B4BFF","#D04939FF","#EB7803FF","#F5A736FF","#FBCF6FFF","#FEF1A6FF","#E2F8B5FF","#9CE5ADFF","#43CBB1FF","#00AAB6FF","#0080B2FF","#584B9FFF")
    }
    vizChord_12x12(data=data,hot=hot,cold=cold,width=width,height=height,filename=filename, colorscheme = colorscheme, colorbar=colorbar)
  } else
  {cat("The length of the input vector does not fit any of the recognized parcellation schemes. The input vector should contain 78, 4005, 7021, 23871 or 30135 values")}
}

########################################################################################################
########################################################################################################
vizChord_edge=function(data,width,height,hot,cold, colorscheme, filename, colorbar)
{

  localenv <- new.env()

  edge_lengths=c(4005,7021,23871,30135)

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

  FC_matrix=array(rep(0,nnode^2),dim=c(nnode,nnode))
  FC_matrix[upper.tri(FC_matrix, diag=FALSE)] = data
  FC_matrix=FC_matrix+t(FC_matrix)

  FCregionmat=array(dim=c(noregions,noregions))
  for (rowno in 1:noregions){
    for (colno in 1:noregions){
      FCregionmat[rowno,colno]=
        mean(FC_matrix
             [label$oldorder[which(regionnoFC==rowno)],label$oldorder[which(regionnoFC==colno)]])
    }
  }

  #region x region matrix FC
  FCregionmat[lower.tri(FCregionmat)]=NA
  x=array(dim=c(noregions^2,3))

  a=1
  for (rowno in 1:noregions){
    for (colno in 1:noregions){
      x[a,1]=rowno
      x[a,2]=colno
      if (rowno==colno) {
        x[a,3]=(FCregionmat[rowno,colno])/2
      }else {
        x[a,3]=FCregionmat[rowno,colno]
      }
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

  if(NROW(which(datFC$value==0))==0){
    datrmFC=datFC
  } else {
    datrmFC=datFC[-which(datFC$value==0),]}

  if(length(unique(c(datrmFC$from,datrmFC$to)))< noregions)
  {
    stop(paste("Connections are absent in",noregions-length(unique(c(datrmFC$from,datrmFC$to))),"of the networks. Please use a connectogram instead",sep=" "))
  }

  for (datval in which(datrmFC$value>0)){
    colarrFC[datval]=pos_color_val[round(abs(datrmFC$value[datval])/max(abs(datrmFC$value))*100)+1]
  }

  for (datval in which(datrmFC$value<0)){
    colarrFC[datval]=neg_color_val[round(abs(datrmFC$value[datval])/max(abs(datrmFC$value))*100)+1]
  }

  FCchord=ggplotify::as.grob(~suppressWarnings(circlize::chordDiagram(datrmFC, col=colarrFC, order=labelnameFC, self.link = 2,
                                                                      grid.col=colorscheme,annotationTrack = c("grid","name"),link.border=colarrFC,)),envir = localenv)
  if(colorbar==T)
  {
    legend.plot=ggplot2::ggplot(datFC, ggplot2::aes(color=value, x=value, y=value))+
      ggplot2::scale_colour_gradient2(name="Connectivity strength",low=cold,mid="white",high=hot,
                                      guide = "colourbar", limits=c(-1,1), breaks=c(-1,1),
                                      labels=c("Strong negative","Strong positive"))+
      ggplot2::geom_point()+
      ggplot2::guides(color=ggplot2::guide_colorbar(ticks=F,title.position = "left",barheight = 0.5))+
      ggplot2::theme(legend.position = "bottom",legend.title=ggplot2::element_text(face="bold", size=7,vjust=1),
                     legend.text =ggplot2::element_text(size=6))
    legend = cowplot::get_legend(legend.plot)

    png(filename =filename, width = width, height = height, res = 300)
    p=cowplot::plot_grid(FCchord,legend,ncol=1, nrow=2, rel_heights = c(0.95,0.05))
    print(p)
  } else if(colorbar==F)
  {
    png(filename =filename, width = width, height = height*0.95, res = 300)
    p=cowplot::plot_grid(FCchord,ncol=1, nrow=1)
    print(p)
  }

  dev.off()
}

########################################################################################################
########################################################################################################

vizChord_12x12=function(data,width,height,hot,cold, colorscheme, filename,colorbar)
{
  localenv = environment()

  nodesname<-c("Auditory","Cinguloopercular","Cinguloparietal",
               "Default\nmode","Dorsal\nattention","Frontoparietal",
               "Retrosplenial\ntemporal","Sensorimotor\nhand","Sensorimotor\nmouth",
               "Salience","Ventral\nattention","Visual")
  no_nodes<-12

  #color parameters
  colarrFC<-array()
  red = colorRampPalette(c("white",hot))
  blue = colorRampPalette(c("white",cold))
  redval<-red(101)
  blueval<-blue(101)

  #reconstruct 12 x 12 FC matrices
  FC_12X12<-array(rep(NA,12^2),dim=c(12,12))
  FC_12X12[upper.tri(FC_12X12, diag=T)] <- data #your results input

  a=1

  datFC0<-array(dim=c(144,3))
  for (rowno in 1:no_nodes){
    for (colno in 1:no_nodes){
      datFC0[a,1]<-rowno
      datFC0[a,2]<-colno
      datFC0[a,3]<-FC_12X12[rowno,colno]
      a=a+1
    }
  }
  #converting 12 x 12 FC matrices to a 'from node to node' data frame
  datFC<-datFC0[-which(is.na(datFC0[,3])),]
  datFC<-data.frame(datFC, stringsAsFactors = F)

  colnames(datFC)<-c("from","to","value")

  #sets up the color function such that darker colors correspond to stronger connections

  col_funFC = circlize::colorRamp2(range(datFC$value), c(hot,cold))

  for (label in 1:no_nodes){
    datFC$from[which(datFC$from==label)]<-nodesname[label]
    datFC$to[which(datFC$to==label)]<-nodesname[label]
  }

  for (datval in which(datFC$value>0)){
    colarrFC[datval]<-redval[round(abs(datFC$value[datval])/max(abs(datFC$value))*100)+1]
  }

  for (datval in which(datFC$value<0)){
    colarrFC[datval]<-blueval[round(abs(datFC$value[datval])/max(abs(datFC$value))*100)+1]
  }

  #generate and saves the FC chord diagram to the FCchord object
  if(colorbar==T)
  {
    FCchord=ggplotify::as.grob(~circlize::chordDiagram(datFC, col=colarrFC, self.link = 1,
                                                       grid.col=colorscheme[1:12], annotationTrack = c("grid","name"),link.border=colarrFC),envir = localenv)
    legend.plot=ggplot2::ggplot(datFC, ggplot2::aes(color=value, x=value, y=value))+
      ggplot2::scale_colour_gradient2(name="Connectivity strength",low=cold,mid="white",high=hot,
                                      guide = "colourbar", limits=c(-1,1), breaks=c(-1,1),
                                      labels=c("Strong negative","Strong positive"))+
      ggplot2::geom_point()+
      ggplot2::guides(color=ggplot2::guide_colorbar(ticks=F,title.position = "left",barheight = 0.5))+
      ggplot2::theme(legend.position = "bottom",legend.title=ggplot2::element_text(face="bold", size=7,vjust=1),
                     legend.text =ggplot2::element_text(size=6))
    legend = cowplot::get_legend(legend.plot)

    png(filename =filename, width = width, height = height, res = 300)
    p=cowplot::plot_grid(FCchord,legend,ncol=1, nrow=2, rel_heights = c(0.95,0.05))
    print(p)
    dev.off()
  } else
  {
    png(filename =filename, width = width, height = height*0.95, res = 300)
    p=cowplot::plot_grid(FCchord,ncol=1, nrow=1)
    print(p)
    dev.off()
  }
}

