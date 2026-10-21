# Generates the neighbour graph for igraph

**\[experimental\]** Symmetrises a kNN matrix into an undirected,
weighted graph. Three options: the bluster/Seurat shared nearest
neighbour graph (`"snn"`), the ScanPy/UMAP fuzzy simplicial set
(`"umap"`) or the plain kNN union (`"knn"`).

## Usage

``` r
rs_sc_knn_graph(knn_mat, knn_dist, graph_params, verbose)
```

## Arguments

- knn_mat:

  Integer matrix. Rows represent cells and the columns represent the
  neighbours (0-indexed!), self excluded and sorted by distance.

- knn_dist:

  Optional numeric matrix of the same shape with the distances. Required
  for `knn_sym = "umap"`.

- graph_params:

  Named list. Recognised keys: `knn_sym` (one of
  `c("snn", "umap", "knn")`); for `"snn"`: `snn_similarity` (one of
  `c("jaccard", "rank")`), `pruning` and `rescue`; for `"umap"`:
  `local_connectivity` and `mix_ratio`; for `"knn"`: `binary`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items:

- edges - Integer vector with the flattened edge pairs (1-indexed!),
  each undirected edge once, ready for
  [`igraph::add_edges()`](https://r.igraph.org/reference/add_edges.html).

- weights - Edge weights of the pairs above.
