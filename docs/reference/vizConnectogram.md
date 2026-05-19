# vizConnectogram

Visualizing brain connectivity profiles with multiple connectogram

## Usage

``` r
vizConnectogram(
  data,
  hot = "#F8766D",
  cold = "#00BFC4",
  ncol = 1,
  nrow = 1,
  edgethickness = 0.8,
  filename = "conn.png",
  colorscheme,
  title,
  width = 1000,
  height = 1050,
  leg.height = 200,
  margins = c(1.2, 1.2, 1, 1.5),
  node.size,
  node.text.size = 1,
  legend.text.size = 7,
  legend.title.size = 8,
  limits,
  title.size = 11,
  title.alignment = "center",
  colorbar_title = "Edge Strength",
  edge_labels = c("Positive", "Negative"),
  row_title
)
```

## Arguments

- data:

  a matrix of edge values with 4005, 7021, 23871 or 30135 columns

- hot:

  color for the positive connections. Set to `#F8766D` by default.

- cold:

  color for the negative connections. Set to `#00BFC4` by default.

- ncol:

  number of columns in the plot. Not used for single row data

- nrow:

  number of rows in the plot. Not used for single row data

- edgethickness:

  a value to adjust the thickness of the edges. Set to 0.8 by default.

- filename:

  output filename with a \*.png file extension. Set to "conn.png" by
  default

- colorscheme:

  an optional vector of color names or color codes to color code the
  networks.

- title:

  a vector of strings to be used as title

- width:

  width (in pixels) of each connectogram. Set to 1000 by default. Not
  used for single row data

- height:

  height (in pixels) of each connectogram . Set to 1050 by default. Not
  used for single row data

- leg.height:

  height (in pixels) of legend, in pixels. Set to 100 by default. Not
  used for single row data

- margins:

  vector of 4 values specifying the amount of empty space on the left,
  right, top and bottom for each connectogram. Set to `c(1.2,1.2,1,1.5)`
  by default. You might want to adjust these values if the text labels
  get cut off by a neighbouring connectogram or the legend. Not used for
  single row data

- node.size:

  size parameter for the dots representing the nodes. If not specified,
  an appropriate size will be set.

- node.text.size:

  size parameter for the text labels on the nodes. Set to 1 by default.

- legend.text.size:

  size parameter for the legend text. Set to 6 by default.

- legend.title.size:

  size parameter for the legend title. Set to 8 by default.

- limits:

  a pair of values that governs the limits of the edge strengths
  displayed. If missing `limits=range(abs(data),na.rm=T)`

- title.size:

  size parameter for the title. Set to 11 by default.

- title.alignment:

  string object specifying title's alignment relatively to its plot
  options are 'left', 'center' (default), 'right'.

- colorbar_title:

  title for the colorbar legend

- edge_labels:

  Vector of two strings defining the labels for the edge legends.
  Default is c("Positive","Negative").

- row_title:

  a vector of strings to be used as left vertical titles for each row of
  plots when there are many

## Value

outputs a .png image

## Details

This function takes a matrix (NROW=number of edges in the connectome;
NCOL=number of edges in the connectome) of edge values and visualizes
the edge-to-edge connectivity with multiple connectograms

## Examples

``` r
if (FALSE) {
data=matrix(sample(c(1,0, -1), 23871*6, replace = TRUE, prob = c(0.001, 0.998,0.001)),nrow=6)
vizConnectogram(data=data,ncol=3, nrow=2)
}
```
