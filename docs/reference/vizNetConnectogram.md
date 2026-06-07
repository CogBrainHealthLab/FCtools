# vizNetConnectogram

Generate a Network Connectogram from a Vector of FC network or edge data

## Usage

``` r
vizNetConnectogram(
  FC_dat,
  show.sig = F,
  title = NULL,
  title.size = 10,
  hot = "#F8766D",
  cold = "#00BFC4",
  node.text.size = 2,
  node.size = 2,
  node.color = "black",
  edge.thickness = 1.5,
  sig.color = "darkgrey",
  legend = T,
  legend.title = "Standardized Coefficient",
  legend.title.size = 6,
  legend.text.size = 5,
  expand = 1.1,
  limits,
  filename,
  width = 1000,
  height = 700
)
```

## Arguments

- show.sig:

  Logical. If `TRUE`, significant edges will be outlined with
  `sig.color`. Default is `FALSE`.

- title:

  Character string for the plot title. Default is `NULL` (no title).

- title.size:

  Numeric. Font size of the plot title. Default is `10`.

- hot:

  Character. Color used for positive edges and node highlights. Default
  is `"#F8766D"` (red).

- cold:

  Character. Color used for negative edges and node highlights. Default
  is `"#00BFC4"` (teal).

- node.text.size:

  Numeric. Size of the node label text. Default is `2`.

- node.size:

  Numeric. Size of the node points. Default is `2`.

- node.color:

  Character. Fill color of the node points. Default is `"black"`.

- edge.thickness:

  Numeric. Width of the edges (and node border stroke). Default is `2`.

- legend:

  Logical. If `TRUE`, a continuous color bar legend is added to the
  plot. Default is `TRUE`.

- legend.title:

  Character string for the legend title. Default is
  `"Standardized Coefficient"`.

- legend.title.size:

  Numeric. Font size of the legend title. Default is `8`.

- legend.text.size:

  Numeric. Font size of the legend tick labels. Default is `5`.

- expand:

  Numeric. Multiplicative expansion factor applied to the x and y axis
  limits to prevent node labels from being clipped. Default is `1.1`.

- limits:

  Numeric vector of length 2. Color scale limits for the edge alpha and
  node alpha mappings. If not supplied, limits are set automatically to
  `c(0, max(abs(weight)))` inside the function.

- filename:

  Character. Path/filename for the output PNG. If the argument is
  omitted entirely, no file is written.

- width:

  Integer. Width of the saved PNG in pixels. Default is `1000`.

- height:

  Integer. Height of the saved PNG in pixels. Default is `1000`.

- data:

  A data frame with row names in the format `"nodeA to nodeB"`, and at
  minimum two columns: `coef` (numeric edge weights / standardized
  coefficients) and `p` (p-values, used when `thresholded = TRUE`).

- show.color:

  characer. Color used for outlining significant edges.Default is
  `darkgrey`.

- thresholded:

  Logical. If `TRUE`, edges with `p > 0.05` are set to zero (i.e.,
  suppressed in the plot). Default is `FALSE`.

## Value

A `ggraph`/`ggplot2` plot object (invisibly returned via `print` and
`return`). As a side effect the plot is printed to the active graphics
device, and — if `filename` is provided — written to a 300 dpi PNG file.

## Details

Creates a circular network connectogram (chord diagram) from a data
frame of connectivity results, visualizing edge weights and their signs
using color and transparency. Node self-connections (diagonal elements)
are highlighted on the node border. Optionally saves the plot to a PNG
file.

Row names of `results` must follow the pattern `"nodeA to nodeB"`.
Self-connections (`from == to`) are extracted to colour and shade node
borders: if all diagonal values are `> 1` the border is drawn in `hot`;
if all are `< 1` it is drawn in `cold`. The diagonal entries are removed
before the `igraph` graph object is constructed.

Edge colour encodes sign (`hot` = positive, `cold` = negative) and edge
transparency encodes absolute magnitude.

The function relies on `igraph`, `ggraph`, and `ggplot2`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example with a dummy results data frame
nodes <- c("A", "B", "C")
pairs <- expand.grid(from = nodes, to = nodes)
rn <- paste(pairs$from, "to", pairs$to)
res <- data.frame(
  coef = rnorm(nrow(pairs)),
  p    = runif(nrow(pairs)),
  row.names = rn
)

p <- vizNetConnectogram(
  results  = res,
  title    = "My Connectogram",
  filename = "connectogram.png"
)
} # }
```
