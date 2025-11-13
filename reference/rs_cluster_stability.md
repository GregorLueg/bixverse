# Helper function to assess cluster stability

Helper function to assess cluster stability

## Usage

``` r
rs_cluster_stability(data)
```

## Arguments

- data:

  Integer matrix. Assumes that each column represents a given
  resampling/bootstrap and the rows represent the features, while each
  integer indicates cluster membership.

## Value

A list containing:

- mean_jaccard - mean Jaccard similarities for this feature across all
  the bootstraps, resamplings.

- std_jaccard - the standard deviation of the Jaccard similarities for
  this feature across all the bootstraps, resamplings.
