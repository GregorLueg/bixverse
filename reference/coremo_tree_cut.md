# Coremo: cuts a hierarchical cluster based on k

This function uses a tree (output of
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)) and cuts it
according to the parameter k. If a `min_size` is specified, modules are
merged by their similarity of their eigen values in the distance matrix.

## Usage

``` r
coremo_tree_cut(
  tree,
  k,
  dist_mat,
  min_size = NULL,
  cor_method = c("pearson", "spearman")
)
```

## Arguments

- tree:

  hclust object. The hierarchical clustering of the correlation matrix
  (or the distance thereof).

- k:

  Integer. Number of cuts on the tree.

- dist_mat:

  Numerical matrix. The distance matrix that was used to compute the
  hierarchical clustering.

- min_size:

  Integer. Optional minimum size for the clusters.

- cor_method:

  String. Which correlation method to use for optionally combining the
  small clusters. One of `c("pearson", "spearman")`.

## Value

A vector with module membership.
