# NBS

Network-based statistics analysis

## Usage

``` r
NBS_lme(
  model,
  contrast,
  random,
  FC_data,
  nperm = 100,
  nthread = 1,
  p = 0.001,
  perm_type = "row"
)
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

- nperm:

  The number of permutations to generate the null distribution of
  network strengths. Set to 100 by default

- nthread:

  The number of CPU threads to use. Set to 1 by default

- p:

  the edge-wise threshold. Set to 0.001 by default

- perm_type:

  A string object specifying whether to permute the rows ("row"),
  between subjects ("between"), within subjects ("within") or between
  and within subjects ("within_between") for random subject effects.
  Default is "row".

## Value

Returns a list object containing

- `results` Edge- and network-wise results in a data.frame object

- `t.orig` Edge-wise t-stats

- `tcrit` The critical t-value

- `max.netstr` A vector containing the null distribution of the permuted
  network strengths

## Details

This function implements the NBS analysis described in [Zalesky et al.
(2010)](https://www.sciencedirect.com/science/article/abs/pii/S1053811910008852))

## Examples

``` r
if (FALSE) { # \dontrun{
model1=NBS_lme(model,contrast, random, FC_data, nperm=1000, nthread=8, p=0.001)
} # }
```
