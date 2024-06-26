% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vizConnectome.R
\name{vizConnectogram}
\alias{vizConnectogram}
\title{vizConnectogram}
\usage{
vizConnectogram(
  data,
  hot = "#F8766D",
  cold = "#00BFC4",
  edgethickness = 0.8,
  filename = "conn.png",
  colorscheme,
  colorbar = T
)
}
\arguments{
\item{data}{a vector of edge values with a length of 4005, 7021, 23871 or 30135.}

\item{hot}{color for the positive connections. Set to \verb{#F8766D} by default.}

\item{cold}{color for the negative connections. Set to \verb{#00BFC4} by default.}

\item{edgethickness}{a value to adjust the thickness of the edges. Set to 0.8 by default.}

\item{filename}{output filename with a *.png file extension. Set to "conn.png" by default}

\item{colorscheme}{an optional vector of color names or color codes to color code the networks.}

\item{colorbar}{\code{TRUE} or \code{FALSE} option on whether a colorbar for the connectivity strength will be displayed. Set to \code{TRUE} by default.}
}
\value{
outputs a .png image
}
\description{
Visualizing brain connectivity profiles with a connectogram
}
\details{
This function takes a vector of edges and visualizes the edge-to-edge connectivity in a connectogram.
}
\examples{

results=sample(c(1,0, -1), 7021, replace = TRUE, prob = c(0.01, 0.98,0.01))
vizConnectogram(data=results, filename="FC_119.png")

}
