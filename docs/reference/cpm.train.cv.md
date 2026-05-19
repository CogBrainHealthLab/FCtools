# cpm.train.cv

K-fold cross-validation procedure to determine the optimal p-value
threshold for feature selection

## Usage

``` r
cpm.train.cv(data, outcome, p, nfolds = 5)
```

## Arguments

- data:

  An N x E matrix containing the vectorized edges; where N = number of
  subjects, E=number of edges

- outcome:

  The outcome variable to predict

- p:

  a vector of p-values to tested. If `p` is not specified, an
  automatically determined geometric sequence of p-values will be tested

- nfolds:

  number of cross-validation folds. Set to 5 by default

## Value

Returns a list object containing

- `opt.pvals` A vector containing the 2 optimal p-values, one each for
  the positive and negative network models

- `results` A P(number of p-values tested) x 2(positive and negative
  network models) matrix containing the predicted-actual correlations
  for each of the p-value selection thresholds in each network model.

- `pvals` A vector of the p-value threshold tested

## Details

This function runs the
[`cpm.train()`](https://cogbrainhealthlab.github.io/FCtools/reference/cpm.train.md)
using a range of p-values and determines the optimal p-values that
maximize predicted-actual correlation in a K-fold cross-validation
paradigm.

## Examples

``` r
if (FALSE) { # \dontrun{
model1.cv=cpm.train.cv(data=FC_data,outcome=dat_beh$age)
} # }
```
