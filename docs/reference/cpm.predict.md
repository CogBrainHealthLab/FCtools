# cpm.predict

predicting scores using a previously trained CPM model

## Usage

``` r
cpm.predict(model, test.data, network = "both")
```

## Arguments

- model:

  A list object generated using
  [`cpm.train()`](https://cogbrainhealthlab.github.io/FCtools/reference/cpm.train.md)

- test.data:

  An N x E matrix containing the vectorized edges of subjects whose
  outcomes/scores are to be predicted. This matrix has to contain the
  same number of columns as the matrix of features used to train the CPM
  model

- network:

  A string to specify which network strength(`'positive'`,`'negative'`
  or `'both'`) to use in the outcome prediction. Set to `'both'` by
  default

## Value

Returns a vector of predicted outcomes if `network='positive'` or
`network='negative'`. In the case of `network='both'`, a data.frame
containing the predicted scores using the positive, negative, and the
combined positive + negative models will be returned

## Details

This function takes a previously trained CPM model's weights to
calculate the positive and negative network strengths for each subject.
Then these network strengths are multiplied with the regression
coefficients in the CPM model to compute a predicted score.

## Examples

``` r
if (FALSE) { # \dontrun{
model1=cpm.predict(model=model1, data=FC_data.test,network="both")
} # }
```
