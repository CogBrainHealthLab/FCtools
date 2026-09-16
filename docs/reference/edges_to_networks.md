# edges_to_networks

convert an edge level FC matrix into a network level FC matrix

## Usage

``` r
edges_to_networks(FCmat)
```

## Arguments

- FCmat:

  an FC matrix or vector

## Value

A network level FC matrix

## Details

This function first identifies the unique network pairs in the
appropriate FC atlas and then averages the edges within each of the
network pairs

## Examples

``` r
edges_to_networks(runif(23871))
#> Error in edges_to_networks(runif(23871)): could not find function "edges_to_networks"
```
