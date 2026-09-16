# intersubject_similarity

This function runs an intersubject similarity analysis to determine if
there is a relationship between similarity in FC and similarity in one
or multiple outcomes.

## Usage

``` r
intersubject_similarity(FC_data, outcome, mode = "diff", nperm = 1000)
```

## Arguments

- FC_data:

  An N x E matrix containing the vectorized edges; where N = number of
  subjects, E=number of edges

- outcome:

  A numerical vector (single outcome) or matrix (multiple outcomes)
  containing the values of the outcome(s) of interest

- mode:

  When set to `"diff"`, FC similarity is calculated as absolute
  difference between the FC vectors of a pair of subjects. When set to
  `"corr"`,FC similarity is calculated as 1 - (pearson's correlation
  coefficient between the FC vectors of a pair of subjects). Set to
  `"diff"` by default.

- nperm:

  number of permutations for the correlation test between FC similarity
  and outcome similarity

## Value

A list object containing

- `FC_difference_matrix` The FC difference matrix

- `Outcome_difference_matrix` The outcome difference matrix

- `permutation_data` The permuted correlation values

## Details

This function runs an intersubject similarity analysis to determine if
there is a relationship between similarity in FC and similarity in one
or multiple outcomes. The outcome(s) will be z-standardized prior to
calculating the intersubject similarity in the outcome(s).

## Examples

``` r
demomat=get('demomat')
results=intersubject_similarity(FC_data = demomat, outcome=c(1,1,2,2),mode="diff")
#> 
#> Correlation between FC and Outcome similarity matrice = 0.328 ; p =0.327
```
