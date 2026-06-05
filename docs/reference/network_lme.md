# network_lme

mass univariate linear mixed effects analysis at the network level

## Usage

``` r
network_lme(model, contrast, random, FC_data, threshold.method = "fdr")
```

## Arguments

- model:

  A data.frame or matrix containing all the predictors in the model

- contrast:

  The predictor of interest. The edge- and network-wise statistics will
  only be estimated for this predictor

- random:

  A N x 1 numeric vector or object containing the values of the random
  variable (optional). Its length should be equal to the number of
  subjects in model (it should NOT be inside the model data.frame).

- FC_data:

  An N x E matrix containing the vectorized edges; where N = number of
  subjects, E=number of edges

- threshold.method:

  method for correcting for multiple tests. set to `fdr` by default

## Value

Returns a data.frame object with `coef` and corrected `p` values

## Details

This function first summarizes the FC edges into their respective
networks and then carry out mass univariate linear mixed effect analyses
on each of the network to network connection

## Examples

``` r
if (FALSE) { # \dontrun{
model1=network_lme(model,contrast, FC_data)
} # }
```
