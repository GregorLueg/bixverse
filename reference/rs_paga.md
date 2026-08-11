# Run PAGA over a kNN graph and a clustering

**\[experimental\]** Implementation of PAGA in Rust. The kNN graph is
read as a directed graph, with an edge running from cell `i` to each of
its neighbours. Distances are not used, hence only the kNN index matrix
is needed here.

## Usage

``` r
rs_paga(knn_mat, partitions, n_partitions)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent cells and the columns represent the
  neighbours (0-indexed!).

- partitions:

  Integer vector. Cluster label per cell (0-indexed!). Needs one entry
  per row of `knn_mat`.

- n_partitions:

  Integer. Declared cluster count. Pass the number of factor levels to
  retain empty ones.

## Value

A list with

- connectivities - The abstracted graph, clusters x clusters, as a
  symmetric CSR list with a zero diagonal. Use
  [`sparse_list_to_mat()`](https://gregorlueg.github.io/bixverse/reference/sparse_list_to_mat.md)
  to transform it into a sparse matrix.

- connectivities_tree - Maximum spanning forest of `connectivities`,
  carrying the original connectivity values on the retained edges. Same
  format as above.

- sizes - Integer vector with the number of cells per cluster.

## References

Wolf, et al., Genome Biol., 2019.
