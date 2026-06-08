# network_lm

mass univariate linear regression at the network level

## Usage

``` r
network_lm(
  model,
  contrast,
  FC_data,
  threshold.method = "fdr",
  perm = T,
  nperm = 1000
)
```

## Arguments

- model:

  A data.frame or matrix containing all the predictors in the model

- contrast:

  The predictor of interest. The edge- and network-wise statistics will
  only be estimated for this predictor

- FC_data:

  An N x E matrix containing the vectorized edges; where N = number of
  subjects, E=number of edges

- threshold.method:

  method for correcting for multiple tests. set to `fdr` by default

- perm:

  If set to `TRUE`, p values will be calculated using a permutation
  approach by shuffling subjects' labels, before correcting for FDR. Set
  to `TRUE` by default

- nperm:

  number of permutations to use if `perm=T`.

## Value

Returns a data.frame object with `coef` and corrected `p` values

## Details

This function first summarizes the FC edges into their respective
networks and then carry out mass univariate linear regression analyses
on each of the network to network connection

## Examples

``` r
if (FALSE) { # \dontrun{
model1=network_lm(model,contrast, FC_data)
} # }
```
