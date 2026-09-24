# Flatten kNN matrix to edge pairs

**\[experimental\]** Helper function to leverage Rust to transform a kNN
matrix into an edge list.

## Usage

``` r
rs_knn_mat_to_edge_pairs(knn_mat, one_index)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent the samples and the columns the 0-based
  indices of the k-nearest neighbours.

- one_index:

  Boolean. Shall 1-based indices be returned.

## Value

A list with the following elements

- from - the from indices

- to - the to indices
