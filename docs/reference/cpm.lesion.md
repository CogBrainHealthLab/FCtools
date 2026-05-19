# cpm.lesion

CPM with the leave-one-network-out lesion approach

## Usage

``` r
cpm.lesion(train.data, test.data, train.outcome, test.outcome, p = 0.05)
```

## Arguments

- train.data:

  An N x E matrix containing the vectorized edges in the training
  dataset; where N = number of subjects, E=number of edges

- test.data:

  An N x E matrix containing the vectorized edges in the training
  dataset

- train.outcome:

  The outcome variable to predict, within the training dataset

- test.outcome:

  The outcome variable to predict, within the testing dataset

- p:

  The p-value threshold of a Pearson's correlation test between the
  feature and outcome that determines if the feature is selected. Set to
  0.05 by default.

## Value

Returns a matrix with the following columns

- `lesion.model` The CPM model in which the listed network is left out

- `positive` The predicted-actual correlations for each of the
  `lesion.model`s' positive network model

- `negative` The predicted-actual correlations for each of the
  `lesion.model`s' negative network model

- `both` The predicted-actual correlations for each of the
  `lesion.model`s' combined positive + negative network model

## Details

This function runs the
[`cpm.train()`](https://cogbrainhealthlab.github.io/FCtools/reference/cpm.train.md)
and
[`cpm.predict()`](https://cogbrainhealthlab.github.io/FCtools/reference/cpm.predict.md)
while leaving out a network of edges each time. The changes in
predicted-actual correlation while a network being left out can be
interpreted to be an indication of the contribution of the network to
predicting the outcome

## Examples

``` r
if (FALSE) { # \dontrun{
model1.lesion=cpm.lesion(train.data=FC_data.train,test.data=FC_data.test,train.outcome=train_dat$age, test.outcome=test_dat$age,p=0.05)
} # }
```
