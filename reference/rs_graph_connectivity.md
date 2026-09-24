# Calculate graph connectivity per cell type

**\[experimental\]** For each cell type, the fraction of its cells in
the largest connected component of the kNN graph restricted to that cell
type. 1 means every cell type forms one connected piece. Edge direction
is ignored.

## Usage

``` r
rs_graph_connectivity(knn_mat, labels)
```

## Arguments

- knn_mat:

  Integer matrix. The rows represent the cells and the columns the
  neighbour indices (0-indexed!).

- labels:

  Integer vector. The cell type per cell. The codes need not be 0-based
  or contiguous.

## Value

A list with the following items

- per_label - Connectivity per cell type, in order of first appearance
  in `labels`.

- mean - Mean connectivity.

- median - Median connectivity.

## References

Luecken, et al., Nat Methods, 2022
