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

a list object containing

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

results=runif(7021, min = -1, max = 1)
vizChord(data=results, filename="FC_chord119.png")


if (FALSE) { # \dontrun{
results=intersub(FC_data = dat_FC, outcome=dat_beh[,10:15],mode="diff")
} # }
```
