# NBS

Network-based statistics analysis

## Usage

``` r
NBS(model, contrast, FC_data, nperm = 100, nthread = 1, p = 0.001)
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

- nperm:

  The number of permutations to generate the null distribution of
  network strengths. Set to 100 by default

- nthread:

  The number of CPU threads to use. Set to 1 by default

- p:

  the edge-wise threshold. Set to 0.001 by default

## Value

A list object containing

- `results` Edge- and network-wise results in a data.frame object

- `t.orig` Edge-wise t-stats

- `tcrit` The critical t-value

- `max.netstr` A vector containing the null distribution of the permuted
  network strengths

## Details

This function implements the NBS analysis described in Zalesky et al.
(2010)
[doi:10.1016/j.neuroimage.2010.06.041](https://doi.org/10.1016/j.neuroimage.2010.06.041)

## Examples

``` r
demomat=get('demomat')[,1:7021] 
contrast=c(1,1,2,2)
random=c('sub1','sub2','sub3','sub4')
model1=NBS(model=contrast, 
           contrast=contrast, 
           FC_data=demomat, 
           nperm=2, 
           nthread=2, 
           p=0.001)
#>   |                                                                              |                                                                      |   0%
#> Estimating permuted network strengths...
#>   |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
#> Completed in :0.1 minutes 
```
