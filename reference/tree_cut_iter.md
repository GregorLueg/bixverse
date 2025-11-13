# Coremo: Iterate over k for gene module detection.

This function uses a tree (output of
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)) and cuts it
according across all values from k_min to k_max and returns the quality
at each individual cut.

## Usage

``` r
tree_cut_iter(
  tree,
  cor_mat,
  dist_mat,
  k_min = 1L,
  k_max = 100L,
  min_size = NULL,
  cor_method = c("spearman", "pearson"),
  seed = 42L
)
```

## Arguments

- tree:

  hclust object. The hierarchical clustering of the correlation matrix
  (or the distance thereof).

- cor_mat:

  Numerical matrix. Correlation matrix.

- dist_mat:

  Numerical matrix. Distance matrix.

- k_min, k_max:

  Integer. The minimum and maximum number of cuts.

- min_size:

  Integer. Optional minimum size of resulting modules.

- cor_method:

  String. Method for the correlation function. One of
  `c("pearson", "spearman")`.

- seed:

  Integer. For reproducibility purposes.

## Value

a data.table with stats (median size of the clusters, median weighted
R^2, and median R^2) on the varying levels of k.
