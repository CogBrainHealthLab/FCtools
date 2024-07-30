## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
#' @title vizHeatmap
#'
#' @description Visualizing brain connectivity matrices with heatmaps
#'
#' @details This function takes a matrix (NROW=number of edges in the connectome; NCOL=number of edges in the connectome) of edge values and visualizes the edge-to-edge connectivity with multiple connectograms
#'
#' @param data a matrix of edge values with 4005, 7021, 23871 or 30135 columns
#' @param ncol number of columns in the plot. Not used for single row data
#' @param nrow number of rows in the plot. Not used for single row data
#' @param hot color for the positive connections. Set to `#F8766D` by default.
#' @param cold color for the negative connections. Set to `#00BFC4` by default.
#' @param colorscheme an optional vector of color names or color codes to color code the networks.
#' @param filename output filename with a *.png file extension. Set to "heatmap.png" by default
#' @param title a vector of strings to be used as title
#' @param width width (in pixels) of each heatmap. Set to 1000 by default. Not used for single row data
#' @param height height (in pixels) of each heatmap . Set to 1000 by default. Not used for single row data
#' @param leg.size height/width (in pixels) of legend, in pixels. Set to 200 by default.
#' @param title.size size parameter for the title. Set to 12 by default.
#' @param legend.title.size size parameter for the legend title. Set to 8 by default.
#' @param legend.text.size size parameter for the legend text. Set to 8 by default.
#' @param line.width line thickness parameter for grid lines dividing the connectome in to networks. Set to 0.3 by default.
#' @param line.color color parameter for grid lines dividing the connectome in to networks. Set to `"black"` by default.
#' @returns outputs a .png image
#'
#' @examples
#' if (FALSE) {
#' data=data=runif(7021,min = -1,max=1)
#' vizHeatmap(data=data)
#' }
#' @importFrom ggplot2 ggplot aes geom_tile scale_y_continuous scale_x_continuous scale_fill_gradient2 ggtitle theme element_line element_text element_blank margin unit
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' @importFrom cowplot plot_grid get_plot_component
#' @importFrom grid grid.raster
#' @importFrom png readPNG
#' @export
########################################################################################################
########################################################################################################
vizHeatmap=function(data,
                    hot="#F8766D", 
                    cold="#00BFC4",
                    ncol,
                    nrow, 
                    filename="heatmap.png", 
                    colorscheme, 
                    line.color="black",
                    title,
                    limits,
                    title.size=12,
                    legend.title.size=8,
                    legend.text.size=7,
                    line.width=0.3,
                    width=1000,
                    height=1000,
                    leg.size=200)
{
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
  nodecol=list(c("#D53E4F","#FC8D59","#FEE08B","#FFFFBF","#E6F598","#99D594","#3288BD"),
                     c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD"))
  nodecol[3:4]=nodecol[2]

  label=labels_dat[[atlas]]
  nodelevels=unique(label$regionlabel)
  label=label[order(label$oldorder),]
  label$regionlabel = factor(label$regionlabel,levels = nodelevels)
  if(missing("colorscheme")){colorscheme = nodecol[[atlas]]}
  if(missing("title")){title=rep("",NROW(data))}
  if(missing("limits")){limits=c(-max(abs(data)),max(abs(data)))}
  limits=c(floor(limits[1]*10)/10,ceiling(limits[2]*10)/10)
  nnodes=NROW(label)
  label=label[order(label$region),]
  
    #determine the end indices of each network
    ends=list()
    count=1
    for(row in 2:NROW(label))
    {
      if(label$regionlabel[row]!=label$regionlabel[row-1])
      {
        ends[[count]]=(row+row-1)/2
        count=count+1
      }
    }
    ends=unlist(ends)
    
    ##dummy data for legend
    dummy_dat=data.frame(1:nnodes,1:nnodes,label$region)
    colnames(dummy_dat)=c("Var1","Var2","value")
    
  ##different steps depending on mode
  if(mode=="vector")
  {
  ##side legend plots
    #node legend
    nodecol=ggplot2::ggplot(data = dummy_dat, ggplot2::aes(x=Var1, y=Var2, fill=as.factor(value))) + 
            ggplot2::geom_tile()+
            ggplot2::scale_fill_manual(values = rev(colorscheme), name="Network",labels=rev(unique(label$regionlabel)))+
            ggplot2::theme(legend.position="right", legend.title.position = "top",
                           legend.text = ggplot2::element_text(size=legend.text.size-2.5),
                           legend.title = ggplot2::element_text(size=legend.title.size-2.5),
                           legend.key.size = ggplot2::unit(0.3,"cm"))
    
    #edge strength legend
    colbar=ggplot2::ggplot(data = dummy_dat, ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
            ggplot2::geom_tile()+
            ggplot2::scale_fill_gradient2(name="Edge strength",limits=round(limits,2),low=cold,high=hot,mid ="white",
                                    na.value = "black",breaks=round(c(limits[1],0,limits[2]),2))+
            ggplot2::theme(legend.position="right",legend.title.position = "top",
                           legend.text = ggplot2::element_text(size=legend.text.size-2.5),
                           legend.title = ggplot2::element_text(size=legend.title.size-2.5),
                           legend.key.size = ggplot2::unit(0.3,"cm"))
  
    #combining both legends          
    legend=cowplot::plot_grid(cowplot::get_plot_component(colbar, 'guide-box-right', return_all = TRUE),
                            cowplot::get_plot_component(nodecol, 'guide-box-right', return_all = TRUE), nrow=2)
    main=genplot_heatmap(row_data=data[1,],title=title[1],nnodes,label, hot, cold,colorscheme,atlas,title.size,limits,ends,line.color,line.width)
    #output plot
    png(filename=filename,width=width+leg.size, height=height, res=300)
    suppressWarnings(print(cowplot::plot_grid(main,legend,ncol=2, rel_widths = c(height,leg.size))))
    dev.off()
    img=png::readPNG(source =filename)
    grid::grid.raster(img)
    
  } else if(mode=="matrix")
  {
    if((nrow*ncol)<NROW(data))  {stop("Not enough columns or rows are specified to display the subplots")}
    #generate multiple plots
    ggplot.obj=list()
    for(row in 1:NROW(data))  {ggplot.obj[[row]]=genplot_heatmap(row_data=data[row,],title=title[row],nnodes,label, hot, cold,colorscheme,atlas,title.size,limits,ends,line.color,line.width)}
    
    ##bottom legend plots
      #node legend
      nodecol=ggplot2::ggplot(data = dummy_dat, ggplot2::aes(x=Var1, y=Var2, fill=as.factor(value))) + 
        ggplot2::geom_tile()+
        ggplot2::scale_fill_manual(values = rev(colorscheme), name="Network",labels=rev(unique(label$regionlabel)))+
        ggplot2::theme(legend.position="bottom", legend.title.position = "top",
                       legend.text = ggplot2::element_text(size=legend.text.size),
                       legend.title = ggplot2::element_text(size=legend.title.size),
                       legend.key.size = ggplot2::unit(0.5,"cm"))
      
      #edge strength legend
      colbar=ggplot2::ggplot(data = dummy_dat, ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
        ggplot2::geom_tile()+
        ggplot2::scale_fill_gradient2(name="Edge strength",limits=round(limits,2),low=cold,high=hot,mid ="white",
                                      na.value = "black",breaks=round(c(limits[1],0,limits[2]),2))+
        ggplot2::theme(legend.position="bottom",legend.title.position = "top",
                       legend.text = ggplot2::element_text(size=legend.text.size),
                       legend.title = ggplot2::element_text(size=legend.title.size),
                       legend.key.size = ggplot2::unit(0.5,"cm"))
    
      #combining both legends
      legend=cowplot::plot_grid(cowplot::get_plot_component(nodecol, 'guide-box-bottom', return_all = TRUE),
                                cowplot::get_plot_component(colbar, 'guide-box-bottom', return_all = TRUE), ncol=2, rel_widths = c(70,30))
      
      #combining legends with heatmaps
      #main=gridExtra::grid.arrange(grobs=ggplot.obj,nrow=nrow,ncol=ncol)
    
      #output plot
      png(filename=filename,width=ncol*width, height=(nrow*height)+leg.size,res=300)
      suppressWarnings(print(cowplot::plot_grid(gridExtra::grid.arrange(grobs=ggplot.obj,nrow=nrow,ncol=ncol),legend,nrow=2, rel_heights = c((nrow*height),leg.size))))
      dev.off()

      img=png::readPNG(source =filename)
      grid::grid.raster(img)
  }
}
  
##function to generate individual plots
genplot_heatmap=function(row_data,title,nnodes,label, hot, cold,colorscheme,atlas,title.size,limits,ends,line.color,line.width)
{
  ##generate first plot
  #reshaping FC vector to FC matrix
  conmat_NxNhalf = matrix(0, nrow = nnodes, ncol = nnodes)
  conmat_NxNhalf[upper.tri(conmat_NxNhalf, diag = FALSE)] = row_data
  conmat_NxN=conmat_NxNhalf+t(conmat_NxNhalf)
  conmat_NxN=conmat_NxN[label$oldorder,label$oldorder]
  
  plotdat=reshape2::melt(conmat_NxN)

  ##graph object
  ggplot2::ggplot(data=plotdat,ggplot2::aes(x=Var1, y=Var2,fill=value))+
    ggplot2::geom_tile()+
    ggplot2::ggtitle(label = title[1])+
    ggplot2::geom_hline(yintercept=ends, color=line.color,linewidth=line.width)+
    ggplot2::geom_vline(xintercept=ends, color=line.color,linewidth=line.width)+
    ggplot2::scale_y_continuous(name="",expand=c(0.005,0), breaks=1:nnodes,labels=label$labels)+
    ggplot2::scale_x_continuous(name="",expand=c(0.005,0), breaks=1:nnodes)+
    ggplot2::scale_fill_gradient2(name="",limits=limits,low=cold,high=hot,mid ="white",
                                  na.value = "black",breaks=c(limits[1],0,limits[2]),guide="none")+
    ggplot2::theme(axis.text=ggplot2::element_blank(),aspect.ratio = 1,
                   axis.ticks = ggplot2::element_line(color=colorscheme[label$region],linewidth = 1,lineend = "square"),
                   plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"), plot.title = ggplot2::element_text(size=title.size))
  
  
}
