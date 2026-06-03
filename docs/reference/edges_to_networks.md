# edges_to_networks

convert an edge level FC matrix into a network level FC matrix

## Usage

``` r
edges_to_networks(FCmat)
```

## Arguments

- FCmat:

  an FC matrix or vector

## Value

returns a network level FC matrix

## Details

This function first identifies the unique network pairs in the
appropriate FC atlas and then averages the edges within each of the
network pairs

## Examples

``` r
edges_to_networks(runif(23871))
#>      Visual to Visual Somatomotor to Visual Somatomotor to Somatomotor
#> [1,]        0.4942361             0.4919609                  0.5019111
#>      Dorsal Attention to Visual Dorsal Attention to Somatomotor
#> [1,]                  0.5016897                       0.4908962
#>      Dorsal Attention to Dorsal Attention Ventral Attention to Visual
#> [1,]                            0.4934034                   0.5148714
#>      Ventral Attention to Somatomotor Ventral Attention to Dorsal Attention
#> [1,]                        0.4943804                             0.5048606
#>      Ventral Attention to Ventral Attention Limbic to Visual
#> [1,]                              0.4993221        0.4939142
#>      Limbic to Somatomotor Limbic to Dorsal Attention
#> [1,]             0.4945089                  0.4911724
#>      Limbic to Ventral Attention Limbic to Limbic Frontoparietal to Visual
#> [1,]                    0.487485        0.4796126                0.4978072
#>      Frontoparietal to Somatomotor Frontoparietal to Dorsal Attention
#> [1,]                     0.4977796                          0.4805301
#>      Frontoparietal to Ventral Attention Frontoparietal to Limbic
#> [1,]                           0.5130455                0.5525792
#>      Frontoparietal to Frontoparietal Default mode to Visual
#> [1,]                        0.5049862               0.506067
#>      Default mode to Somatomotor Default mode to Dorsal Attention
#> [1,]                   0.5067924                        0.4873407
#>      Default mode to Ventral Attention Default mode to Limbic
#> [1,]                         0.4874604              0.4723027
#>      Default mode to Frontoparietal Default mode to Default mode
#> [1,]                       0.498743                    0.5066932
#>      Subcortical to Visual Subcortical to Somatomotor
#> [1,]             0.5062473                  0.4863675
#>      Subcortical to Dorsal Attention Subcortical to Ventral Attention
#> [1,]                           0.473                        0.4909915
#>      Subcortical to Limbic Subcortical to Frontoparietal
#> [1,]             0.4730012                     0.4797091
#>      Subcortical to Default mode Subcortical to Subcortical
#> [1,]                    0.512936                  0.5031385
```
