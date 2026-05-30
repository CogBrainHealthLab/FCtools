#' @title vizNetConnectogram
#' @description Generate a Network Connectogram from a Vector of Results
#' 
#'
#' @details Creates a circular network connectogram (chord diagram) from a data frame of
#' connectivity results, visualizing edge weights and their signs using color and
#' transparency. Node self-connections (diagonal elements) are highlighted on the
#' node border. Optionally saves the plot to a PNG file.
#' @param results A data frame with row names in the format \code{"nodeA to nodeB"},
#'   and at minimum two columns: \code{coef} (numeric edge weights / standardized
#'   coefficients) and \code{p} (p-values, used when \code{thresholded = TRUE}).
#' @param thresholded Logical. If \code{TRUE}, edges with \code{p > 0.05} are set
#'   to zero (i.e., suppressed in the plot). Default is \code{FALSE}.
#' @param title Character string for the plot title. Default is \code{NULL} (no title).
#' @param title.size Numeric. Font size of the plot title. Default is \code{10}.
#' @param hot Character. Color used for positive edges and node highlights.
#'   Default is \code{"#F8766D"} (red).
#' @param cold Character. Color used for negative edges and node highlights.
#'   Default is \code{"#00BFC4"} (teal).
#' @param node.text.size Numeric. Size of the node label text. Default is \code{2}.
#' @param node.size Numeric. Size of the node points. Default is \code{2}.
#' @param node.color Character. Fill color of the node points. Default is
#'   \code{"black"}.
#' @param edge.thickness Numeric. Width of the edges (and node border stroke).
#'   Default is \code{2}.
#' @param legend Logical. If \code{TRUE}, a continuous color bar legend is added
#'   to the plot. Default is \code{TRUE}.
#' @param legend.title Character string for the legend title. Default is
#'   \code{"Standardized Coefficient"}.
#' @param legend.title.size Numeric. Font size of the legend title. Default is
#'   \code{8}.
#' @param legend.text.size Numeric. Font size of the legend tick labels. Default
#'   is \code{5}.
#' @param expand Numeric. Multiplicative expansion factor applied to the x and y
#'   axis limits to prevent node labels from being clipped. Default is \code{1.1}.
#' @param limits Numeric vector of length 2. Color scale limits for the edge alpha
#'   and node alpha mappings. If not supplied, limits are set automatically to
#'   \code{c(0, max(abs(weight)))} inside the function.
#' @param filename Character. Path/filename for the output PNG. If the argument is
#'   omitted entirely, no file is written.
#' @param width Integer. Width of the saved PNG in pixels. Default is \code{1000}.
#' @param height Integer. Height of the saved PNG in pixels. Default is \code{1000}.
#'
#' @return A \code{ggraph}/\code{ggplot2} plot object (invisibly returned via
#'   \code{print} and \code{return}). As a side effect the plot is printed to the
#'   active graphics device, and — if \code{filename} is provided — written to a
#'   300 dpi PNG file.
#'
#' @details
#' Row names of \code{results} must follow the pattern \code{"nodeA to nodeB"}.
#' Self-connections (\code{from == to}) are extracted to colour and shade node
#' borders: if all diagonal values are \code{> 1} the border is drawn in
#' \code{hot}; if all are \code{< 1} it is drawn in \code{cold}. The diagonal
#' entries are removed before the \code{igraph} graph object is constructed.
#'
#' Edge colour encodes sign (\code{hot} = positive, \code{cold} = negative) and
#' edge transparency encodes absolute magnitude.
#'
#' The function relies on \code{igraph}, \code{ggraph}, and \code{ggplot2}.
#'
#' @examples
#' \dontrun{
#' # Minimal example with a dummy results data frame
#' nodes <- c("A", "B", "C")
#' pairs <- expand.grid(from = nodes, to = nodes)
#' rn <- paste(pairs$from, "to", pairs$to)
#' res <- data.frame(
#'   coef = rnorm(nrow(pairs)),
#'   p    = runif(nrow(pairs)),
#'   row.names = rn
#' )
#'
#' p <- vizNetConnectogram(
#'   results  = res,
#'   title    = "My Connectogram",
#'   filename = "connectogram.png"
#' )
#' }
#'
#' @importFrom igraph graph_from_data_frame V E edge_attr edge_attr<-
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text
#'   expand_limits theme_graph
#' @importFrom ggplot2 aes ggtitle geom_point scale_fill_gradientn
#'   scale_edge_color_manual scale_edge_alpha_continuous scale_alpha_continuous
#'   guides guide_colorbar theme element_text unit
#' @export
############################################################################################################################
############################################################################################################################

vizNetConnectogram=function(results,
                                      thresholded=F,
                                      title=NULL,
                                      title.size=10,
                                      hot="#F8766D", 
                                      cold="#00BFC4",
                                      node.text.size=2, 
                                      node.size=2,
                                      node.color="black",
                                      edge.thickness=2,
                                      legend=T,
                                      legend.title="Standardized Coefficient",
                                      legend.title.size=8,
                                      legend.text.size=5,
                                      expand=1.1,
                                      limits,
                                      filename,
                                      width=1000, 
                                      height=1000)
{
  weight=results$coef
  if(thresholded==T)  {weight[results$p>0.05]=0}
  
  #generate graph object
  graph.dat=data.frame(from=sub(" to .*", "", rownames(results)),
                       to=sub(".* to ", "", rownames(results)),
                       weight=weight)
  
  diag=graph.dat$weight[which(graph.dat$from==graph.dat$to)]
  diag=diag[c(2:8,1)] ##need to reoder diagonal values
  diag.alpha=diag
  diag.color=diag
  if(all(diag>1))
  {
    diag.color=hot
  } else if(all(diag<1))
  {
    diag.color=cold
  }
  
  #set graph parameter
  graph.dat=graph.dat[-which(graph.dat$from==graph.dat$to),]
  graph.obj=igraph::graph_from_data_frame(graph.dat, directed = F)
  
  n_nodes=length(igraph::V(graph.obj)$name)
  
  edge.weight=igraph::edge_attr( graph.obj, "weight", index = igraph::E(graph.obj))
  edge.sign=edge.weight
  
  edge.sign=replace(edge.sign, which(edge.sign < 0), "neg")
  edge.sign=replace(edge.sign, which(edge.sign!="neg"), "pos")
  
  igraph::edge_attr(graph.obj, "weight", index = igraph::E(graph.obj))=abs(as.numeric(edge.weight))
  igraph::edge_attr(graph.obj, "signposneg", index = igraph::E(graph.obj))=edge.sign
  limits=c(0,max(abs(weight)))
  
  #plot
  plot.obj=ggraph(graph.obj,layout = 'linear', circular = TRUE)+
    ggtitle(title)+
    geom_node_text(ggplot2::aes(label=igraph::V(graph.obj)$name, x = x * 1.1, y = y* 1.1,
                                angle = ifelse(atan(-(x/y))*(180/pi) < 0,atan(-(x/y))*(180/pi), atan(-x/y)*(180/pi)),
                                hjust = 0.5),size=node.text.size) +
    expand_limits(x =c(-expand,expand), y = c(-expand,expand))+
    geom_edge_arc(aes(color=signposneg,alpha=weight),width=edge.thickness,show.legend = F)+
    geom_node_point(color=node.color, size=node.size,show.legend = F)+
    geom_node_point(shape = 21,aes(color=rep(diag.color,8), alpha=abs(diag.alpha)), size=5,stroke=edge.thickness, show.legend = F)+
    scale_edge_color_manual(guide="none",values = c("pos"=hot,"neg"=cold))+
    scale_edge_alpha_continuous(guide="none",range=c(0,1),limits=limits)+
    scale_alpha_continuous(guide="none",range=c(0,1), limits=limits)+
    theme_graph(background = 'white', text_colour = 'black', bg_text_colour = 'black')+
    theme(plot.title = element_text(size=title.size),aspect.ratio = 1, plot.margin=unit(c(0,0,1,0),units = "cm"))
  
  if(legend==T)
  {
    plot.obj=plot.obj+
      geom_point(data = data.frame(x = NA_real_, y = NA_real_, fill_val = NA_real_),aes(x = x, y = y, fill = fill_val),na.rm = TRUE, inherit.aes = FALSE) +
      scale_fill_gradientn(colours = c(cold,"white",hot),limits= c(-abs(limits[2]),abs(limits[2])), name= legend.title, na.value = NA) +
      guides(fill = guide_colorbar(title.position = "left",
                                   barwidth  = 0.75, 
                                   theme = theme(legend.text = element_text(size = legend.text.size)),
                                   title.theme = element_text(angle = 90, size=legend.title.size)))
  }
  
  if(!missing(filename))
  {
    png(filename=filename, width=width, height=height, res=300)
    print(plot.obj)
    dev.off()
  }
  print(plot.obj)
  return(plot.obj)
}

