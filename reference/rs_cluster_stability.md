# Helper function to assess cluster stability

**\[experimental\]**

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

A list containing, one entry per feature:

- mean_jaccard - mean Jaccard similarity of the feature's cluster across
  all pairs of bootstraps/resamplings.

- std_jaccard - the (population) standard deviation of these Jaccard
  similarities.
