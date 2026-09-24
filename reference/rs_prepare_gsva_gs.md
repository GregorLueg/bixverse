# Prepare a pathway list for GSVA

**\[experimental\]**

## Usage

``` r
rs_prepare_gsva_gs(feature_names, pathway_list, min_size, max_size)
```

## Arguments

- feature_names:

  String vector. The feature names of the matrix (should be the rows).

- pathway_list:

  List. A list containing the pathways and the respective genes

- min_size, max_size:

  Integer. The minimum and maximum size respectively.

## Value

Returns a named list with the sorted, 0-based indices of the pathways
whose size (after matching to `feature_names`) is within
`[min_size, max_size]`.
