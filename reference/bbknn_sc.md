# Run BBKNN

This function implements the batch-balanced k-nearest neighbour
algorithm from Polański, et al. Briefly, the algorithm generate a KNN
index on a per batch basis and identifies the neighbours of cells for
each individual index. Subsequently, it leverages UMAP connectivity
calculations to reduce spurious connections. For more details, please
refer to Polański, et al.

## Usage

``` r
bbknn_sc(
  object,
  batch_column,
  no_neighbours_to_keep = 15L,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  bbknn_params = params_sc_bbknn(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- no_neighbours_to_keep:

  Integer. Maximum number of neighbours to keep from the BBKNN
  algorithm. Due to generating neighbours for each batch, there might be
  a large number of generated neighbours. This will only keep the top
  `no_neighbours_to_keep` neighbours.

- embd_to_use:

  String. The embedding to use. Atm, the only option is `"pca"`.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- bbknn_params:

  A list, please see [`params_sc_bbknn()`](params_sc_bbknn.md). The list
  has the following parameters:

  - neighbours_within_batch - Integer. Number of neighbours to consider
    per batch.

  - knn_method - String. One of `c("annoy", "hnsw")`. Defaults to
    `"annoy"`.

  - ann_dist - String. One of `c("cosine", "euclidean")`. The distance
    metric to be used for the approximate neighbour search. Defaults to
    `"cosine"`.

  - set_op_mix_ratio - Numeric. Mixing ratio between union (1.0) and
    intersection (0.0).

  - local_connectivity - Numeric. UMAP connectivity computation
    parameter, how many nearest neighbours of each cell are assumed to
    be fully connected.

  - annoy_n_trees - Integer. Number of trees to use in the generation of
    the Annoy index.

  - search_budget - Integer. Search budget per tree for the `annoy`
    algorithm.

  - trim - Optional integer. Trim the neighbours of each cell to these
    many top connectivities. May help with population independence and
    improve the tidiness of clustering. If `NULL`, it defaults to
    `10 * neighbours_within_batch`.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The object with added kNN matrix based on BBKNN and the graph based on
the returned connectivities of the algorithm.

## References

Polański, et al., Bioinformatics, 2020
