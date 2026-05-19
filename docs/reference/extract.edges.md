# extract.edges

Generates edge-wise masks for calculating subject-level network
strengths

## Usage

``` r
extract.edges(NBS.obj, network = 1)
```

## Arguments

- NBS.obj:

  A list object generated from an earlier
  [`NBS()`](https://cogbrainhealthlab.github.io/FCtools/reference/NBS.md)
  analysis

- network:

  the network number (reported in the earlier NBS results) of the
  network to be masked. Set to 1 by default

## Value

Returns a list object containing

- `clust.tstat` thresholded edge-wise t-statistics.Edges not belonging
  to this cluster will be zeroed.

- `pos.edges` A vector of 1s and 0s indicating the significant
  network-thresholded positive edges.

- `neg.edges` A vector of -1s and 0s indicating the significant
  network-thresholded negative edges.

- `pos.mask` A vector of 1s and 0s indicating the significant
  network-thresholded positive edges.

- `neg.mask` A vector of 1s and 0s indicating the significant
  network-thresholded negative edges.

## Details

This function generates positive and negative masks (vectors of 1s and
0s), where 1s indicate a significant network-thresholded edge. These
masks can then be used to perform a matrix multiplication with the
vectorized FC matrices to object subject-level network strengths

## Examples

``` r
if (FALSE) { # \dontrun{
extract.edges(model1,network=1)
} # }
```
