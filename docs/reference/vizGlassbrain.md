# vizGlassbrain

Visualizing brain connectivity with a glass brain plot

## Usage

``` r
vizGlassbrain(
  data,
  surf_color = "grey",
  node_color = "#00BA38",
  surf_alpha = 0.2,
  cmap,
  node_size = 8,
  node_label = TRUE,
  node_label_size = 10,
  node_label_color = "black",
  edgethickness = 8,
  limits,
  colorbar_title = "Connectivity strength",
  orientation_labels = TRUE,
  remove_brain = FALSE
)
```

## Arguments

- data:

  a vector of edge values with a length of 4005, 7021, 23871 or 30135

- surf_color:

  color of the cortical surface. Set to `'grey'` by default

- node_color:

  color of the nodes. Set to `'#00BA38'` by default

- surf_alpha:

  alpha value of the cortical surface, where 0 will cause the cortical
  surface to disappear and 1 will cause the cortical surface to be
  completely opaque. Set to `0.2` by default

- cmap:

  A string vector containing 2 to 3 color names/codes specifying the
  colors to be used for the color scale. See
  [`RColorBrewer::display.brewer.all()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  for all possible cmap options. If none are specified, appropriate
  colors will be automatically selected according to `range(data)`

- node_size:

  size parameter for the dots representing the nodes. Set to `8` by
  default.

- node_label:

  option to show node labels. Set to `TRUE` by default.

- node_label_size:

  font size of the node labels. Set to `10` by default.

- node_label_color:

  font color of the node labels. Set to `black` by default.

- edgethickness:

  a value to adjust the thickness of the edges. Set to `8` by default.

- limits:

  Numeric vector of length 2, composed of the lower and upper color
  scale limits of the plot. When left unspecified, the symmetrical
  limits `c(-max(abs(data),max(abs(data)))` will be used.

- colorbar_title:

  title for the colorbar legend. Set to `'Connectivity strength` by
  default

- orientation_labels:

  A boolean object specifying if orientation labels are to be displayed.
  Set to `TRUE` by default

- remove_brain:

  A boolean object specifying cortical surface should be removed. Set to
  `FALSE` by default

## Value

A plot_ly object

## Details

This function takes a vector of edge values and visualizes the
edge-to-edge connectivity with a glass brain plot

## Examples

``` r
mask=sample(c(1,0), 7021, replace = TRUE, prob = c(0.001, 0.999))
data=runif(7021,min = -1,max=1)*mask
vizGlassbrain(data,orientation_labels = TRUE)
#> Error in igraph::graph_from_adjacency_matrix(FCmat.bin): Cannot create a graph object because the adjacency matrix contains NAs.
```
