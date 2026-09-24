# Flatten kNN matrix to edge list

**\[experimental\]** Helper function to leverage Rust to transform a kNN
matrix into a flat edge list.

## Usage

``` r
rs_knn_mat_to_edge_list(knn_mat, one_index)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent the samples and the columns the 0-based
  indices of the k-nearest neighbours.

- one_index:

  Boolean. Shall 1-based indices be returned.

## Value

A flat vector representing the edge list, alternating from and to, i.e.
`c(from_1, to_1, from_2, to_2, ...)`.
