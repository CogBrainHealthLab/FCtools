# network_lm

mass univariate linear regression at the network level

## Usage

``` r
network_lm(model, contrast, FC_data, threshold.method = "fdr")
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

## Value

Returns a data.frame object with `coef` and corrected `p` values

## Details

This function implements the NBS analysis described in [Zalesky et al.
(2010)](https://www.sciencedirect.com/science/article/abs/pii/S1053811910008852))

## Examples

``` r
if (FALSE) { # \dontrun{
model1=network_lm(model,contrast, FC_data)
} # }
```
