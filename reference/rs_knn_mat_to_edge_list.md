# Flatten kNN matrix to edge list

Helper function to leverage Rust to transform a kNN matrix into an edge
list.

## Usage

``` r
rs_knn_mat_to_edge_list(knn_mat, one_index)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent the samples and the columns the indices
  of the k-nearest neighbours.

- one_index:

  Boolean. If the original data is 0-index, shall 1-indexed data be
  returned.

## Value

A flat vector representing the edge list.
