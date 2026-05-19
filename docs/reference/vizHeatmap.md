# vizHeatmap

Visualizing brain connectivity matrices with heatmaps

## Usage

``` r
vizHeatmap(
  data,
  hot = "#F8766D",
  cold = "#00BFC4",
  mid = "white",
  ncol,
  nrow,
  filename = "heatmap.png",
  colorscheme,
  line.color = "black",
  title,
  limits,
  title.size = 12,
  legend.title.size = 8,
  legend.text.size = 7,
  line.width = 0.3,
  width = 1000,
  height = 1000,
  leg.size = 200,
  colorbar_title = "Edge Strengt"
)
```

## Arguments

- data:

  a matrix of edge values with 4005, 7021, 23871 or 30135 columns

- hot:

  color for the positive connections. Set to `#F8766D` by default.

- cold:

  color for the negative connections. Set to `#00BFC4` by default.

- mid:

  color for the mid point of the color scale. Set to `white` by default.

- ncol:

  number of columns in the plot. Not used for single row data

- nrow:

  number of rows in the plot. Not used for single row data

- filename:

  output filename with a \*.png file extension. Set to "heatmap.png" by
  default

- colorscheme:

  an optional vector of color names or color codes to color code the
  networks.

- line.color:

  color parameter for grid lines dividing the connectome in to networks.
  Set to `"black"` by default.

- title:

  a vector of strings to be used as title

- title.size:

  size parameter for the title. Set to 12 by default.

- legend.title.size:

  size parameter for the legend title. Set to 8 by default.

- legend.text.size:

  size parameter for the legend text. Set to 8 by default.

- line.width:

  line thickness parameter for grid lines dividing the connectome in to
  networks. Set to 0.3 by default.

- width:

  width (in pixels) of each heatmap. Set to 1000 by default. Not used
  for single row data

- height:

  height (in pixels) of each heatmap . Set to 1000 by default. Not used
  for single row data

- leg.size:

  height/width (in pixels) of legend, in pixels. Set to 200 by default.

- colorbar_title:

  title for the colorbar legend

## Value

outputs a .png image

## Details

This function takes a matrix (NROW=number of edges in the connectome;
NCOL=number of edges in the connectome) of edge values and visualizes
the edge-to-edge connectivity with multiple connectograms

## Examples

``` r
if (FALSE) {
data=data=runif(7021,min = -1,max=1)
vizHeatmap(data=data)
}
```
