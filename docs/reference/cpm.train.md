# cpm.train

Training a connectome-based predictive model using a fixed p-value
threshold

## Usage

``` r
cpm.train(data, outcome, p = 0.05)
```

## Arguments

- data:

  An N x E matrix containing the vectorized edges; where N = number of
  subjects, E=number of edges

- outcome:

  The outcome variable to predict

- p:

  The p-value threshold of a Pearson's correlation test between the
  feature and outcome that determines if the feature is selected. Set to
  0.05 by default.

## Value

Returns a list object containing

- `weights` vector of 1s (positive edges), 0s and -1s (negative edges)
  indicating the selected edges

- `pos.network.coef` linear regression coefficients of the positive
  network strength model

- `neg.network.coef` linear regression coefficients of the negative
  network strength model

- `both.network.coef` linear regression coefficients of the positive and
  negative network strength model

## Details

This function implements the connectome-based prediction model described
in [Shen et al. (2017)](https://www.nature.com/articles/nprot.2016.178))

## Examples

``` r
if (FALSE) { # \dontrun{
model1=cpm.train(data=FC_data,outcome=dat_beh$age, p=0.05)
} # }
```
