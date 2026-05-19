# edgelist

Listing the edges in a connectogram.

## Usage

``` r
edgelist(data)
```

## Arguments

- data:

  a vector of edge values with a length of 4005, 7021, 23871 or 30135.

## Value

a data.frame() object with the columns of `node_1`, `node_2` and
`weight` columns

## Details

This function takes a vector of edge-to-edge connections, and returns a
data.frame object containing the edges (defined by their pair of
connecting nodes) and their edge values

## Examples

``` r

results=sample(c(1,0, -1), 7021, replace = TRUE, prob = c(0.01, 0.98,0.01))
edgelist(data=results)
#>                     node_1                 node_2 weight
#> 1                  R_Vis_1          R_FrOperIns_1      1
#> 2                  R_Vis_2             R_SomMot_1     -1
#> 3                  R_Vis_2                R_Med_2      1
#> 4                  R_Vis_3               L_Cing_1     -1
#> 5                  R_Vis_3               L_PFCl_1     -1
#> 6                  R_Vis_3                L_Vis_3      1
#> 7                  R_Vis_4           L_TempPole_2     -1
#> 8                  R_Vis_4               L_Post_3      1
#> 9                  R_Vis_5                R_OFC_1      1
#> 10                 R_Vis_7               R_PFCv_1     -1
#> 11                 R_Vis_7              R_Putamen     -1
#> 12                 R_Vis_7                L_PFC_2     -1
#> 13                 R_Vis_7                L_Vis_4     -1
#> 14                 R_Vis_8               L_PFCl_1      1
#> 15              R_SomMot_1               R_PFCl_4     -1
#> 16              R_SomMot_1              L_Caudate     -1
#> 17              R_SomMot_1               L_PFCl_1     -1
#> 18              R_SomMot_1             L_SomMot_5      1
#> 19              R_SomMot_2              L_Putamen     -1
#> 20              R_SomMot_3            _Brain_Stem      1
#> 21              R_SomMot_4                R_Par_2     -1
#> 22              R_SomMot_4                L_Vis_2     -1
#> 23              R_SomMot_5                L_PFC_4      1
#> 24              R_SomMot_5            L_ParOper_1     -1
#> 25              R_SomMot_5             L_SomMot_3     -1
#> 26              R_SomMot_5                L_Vis_5     -1
#> 27              R_SomMot_6           L_Cerebellum      1
#> 28              R_SomMot_7               L_Post_5      1
#> 29                R_Post_1                L_Med_1      1
#> 30                R_Post_1                L_Vis_5      1
#> 31                R_Post_2               R_PFCl_2      1
#> 32                R_Post_2                L_PFC_2      1
#> 33                R_Post_2               L_PrCv_1      1
#> 34                R_Post_3               R_PFCl_3      1
#> 35                R_Post_3               L_Temp_2     -1
#> 36                R_Post_4                L_PFC_6      1
#> 37                R_Post_4                L_PFC_5      1
#> 38                 R_FEF_1              R_PFCmp_1     -1
#> 39          R_TempOccPar_2            R_pCunPCC_2     -1
#> 40          R_TempOccPar_2               L_Post_3      1
#> 41          R_TempOccPar_2               L_Post_1      1
#> 42           R_FrOperIns_1               R_PFCv_1     -1
#> 43           R_FrOperIns_1          L_Hippocampus     -1
#> 44           R_FrOperIns_1               L_PFCl_1     -1
#> 45                 R_Med_1                R_Par_2     -1
#> 46                 R_Med_1                L_Med_3     -1
#> 47                 R_Med_1             L_SomMot_5     -1
#> 48            R_TempPole_1             L_SomMot_2      1
#> 49                 R_Par_1          L_FrOperIns_2      1
#> 50                 R_Par_2            _Brain_Stem     -1
#> 51                R_PFCl_1 L_Diencephalon_Ventral     -1
#> 52                R_PFCl_2               R_PFCv_1      1
#> 53                R_PFCl_2            L_pCunPCC_1     -1
#> 54                R_PFCl_2                L_Vis_6      1
#> 55                R_PFCl_3 L_Diencephalon_Ventral     -1
#> 56                R_PFCl_3                L_Vis_7     -1
#> 57                R_PFCl_4               R_PFCv_1      1
#> 58                R_Cing_1                L_Vis_1     -1
#> 59               R_PFCmp_1               R_PFCv_2     -1
#> 60               R_PFCmp_1             L_Thalamus     -1
#> 61               R_PFCmp_1                L_Vis_8     -1
#> 62                R_pCun_1            R_pCunPCC_1     -1
#> 63                R_pCun_1            L_ParOper_1      1
#> 64                 R_Par_1          R_Hippocampus      1
#> 65                 R_Par_1            L_pCunPCC_2     -1
#> 66                 R_Par_1                L_Vis_9      1
#> 67                R_Temp_1             R_Amygdala     -1
#> 68                R_Temp_2           R_PFCdPFCm_1      1
#> 69                R_Temp_2             L_SomMot_3      1
#> 70                R_PFCv_1            R_pCunPCC_2     -1
#> 71                R_PFCv_1                L_PFC_1     -1
#> 72                R_PFCv_1            L_ParOper_1     -1
#> 73                R_PFCv_1               L_Post_2      1
#> 74                R_PFCv_2            R_pCunPCC_2      1
#> 75                R_PFCv_2             L_SomMot_6      1
#> 76            R_PFCdPFCm_1               L_Temp_2     -1
#> 77            R_PFCdPFCm_1          L_FrOperIns_2     -1
#> 78            R_PFCdPFCm_2           L_TempPole_1     -1
#> 79            R_PFCdPFCm_3              L_Putamen     -1
#> 80            R_PFCdPFCm_3                L_Med_3     -1
#> 81             R_pCunPCC_2             L_Amygdala      1
#> 82             R_Accumbens             R_Thalamus      1
#> 83             R_Accumbens             L_Pallidum      1
#> 84             R_Accumbens L_Diencephalon_Ventral     -1
#> 85             R_Accumbens                L_PFC_1      1
#> 86              R_Amygdala          L_Hippocampus     -1
#> 87           R_Hippocampus              R_Putamen      1
#> 88           R_Hippocampus               L_Post_6     -1
#> 89              R_Pallidum          L_Hippocampus     -1
#> 90               R_Putamen             L_Pallidum     -1
#> 91               R_Putamen L_Diencephalon_Ventral      1
#> 92              R_Thalamus            L_pCunPCC_2      1
#> 93            R_Cerebellum                L_Vis_6     -1
#> 94             _Brain_Stem                L_Vis_4      1
#> 95            L_Cerebellum                L_PFC_3     -1
#> 96            L_Cerebellum               L_Cing_1      1
#> 97            L_Cerebellum                L_Vis_7      1
#> 98              L_Thalamus             L_SomMot_4      1
#> 99               L_Putamen             L_SomMot_3     -1
#> 100             L_Pallidum          L_FrOperIns_1     -1
#> 101             L_Amygdala               L_PFCl_1     -1
#> 102             L_Amygdala                L_Vis_8      1
#> 103             L_Amygdala                L_Vis_5      1
#> 104            L_Accumbens           L_TempPole_2     -1
#> 105            L_Accumbens                L_Med_1      1
#> 106 L_Diencephalon_Ventral                L_PFC_4     -1
#> 107                L_PFC_6                L_Vis_7      1
#> 108                L_PFC_6                L_Vis_1      1
#> 109                L_PFC_5                L_FEF_1      1
#> 110                L_PFC_4                L_Par_1      1
#> 111                L_PFC_2           L_TempPole_2     -1
#> 112                L_Par_2          L_FrOperIns_1     -1
#> 113               L_Temp_1                L_Vis_7     -1
#> 114               L_Cing_1                L_Vis_9     -1
#> 115               L_PFCl_1                L_Med_1     -1
#> 116                L_Par_1             L_SomMot_6     -1
#> 117                L_Par_1                L_Vis_8     -1
#> 118           L_TempPole_1                L_Med_2      1
#> 119                L_Med_3                L_Vis_3      1
#> 120                L_Med_1             L_SomMot_5      1
#> 121               L_PFCl_1            L_ParOper_1     -1
#> 122               L_PrCv_1             L_SomMot_2     -1
#> 123               L_Post_6               L_Post_1     -1
#> 124               L_Post_5               L_Post_4      1
#> 125               L_Post_4                L_Vis_3      1
#> 126               L_Post_2             L_SomMot_3      1
#> 127             L_SomMot_6             L_SomMot_5      1
#> 128             L_SomMot_5                L_Vis_3      1
#> 129             L_SomMot_1                L_Vis_2      1
#> 130                L_Vis_6                L_Vis_5     -1
```
