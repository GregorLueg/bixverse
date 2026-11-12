# Generate a weighted nearest neighbour (WNN) graph

This function implements the approach from Hao et al., to generate a
weighted nearest neighbour graph given two 'omics modalities (for
example RNA and ADT). Per-cell modality weights are computed by
comparing within- and cross-modality neighbourhood distances, and these
weights are used to fuse the two modalities into a single multimodal kNN
graph.

The per-modality kNN graphs are always recomputed internally from the
chosen embeddings at `knn_range` neighbours (see
[`params_sc_wnn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_wnn.md));
the kNN graphs previously stored on the object via
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md)
are not reused, as WNN requires a larger candidate pool than the default
neighbour search.

The result is stored in the object's `other_data` under `"wnn"`,
containing the fused kNN graph as a `SingleCellNearestNeighbour` (with a
kernelised pseudo-distance metric), the neighbour graph as an igraph and
a table of per-cell modality weights. This graph can subsequently be
used for clustering and 2D embedding via the `"wnn"` modality.

## Usage

``` r
generate_wnn_graph_sc(
  object,
  modality_1 = "rna",
  modality_2 = "adt",
  embd_to_use_1 = "pca",
  embd_to_use_2 = "pca",
  no_embd_to_use_1 = NULL,
  no_embd_to_use_2 = NULL,
  wnn_params = params_sc_wnn(),
  graph_params = params_knn_graph_defaults(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class to which to add the WNN.

- modality_1:

  String. First modality. Defaults to `"rna"`.

- modality_2:

  String. Second modality. Defaults to `"adt"`.

- embd_to_use_1:

  String. The embedding to use for the first modality. Must be available
  in the object for `modality_1`. Defaults to `"pca"`.

- embd_to_use_2:

  String. The embedding to use for the second modality. Must be
  available in the object for `modality_2`. Defaults to `"pca"`.

- no_embd_to_use_1:

  Optional integer. Number of embedding dimensions to use for
  `embd_to_use_1`. If `NULL`, all will be used.

- no_embd_to_use_2:

  Optional integer. Number of embedding dimensions to use for
  `embd_to_use_2`. If `NULL`, all will be used.

- wnn_params:

  Named list. Controls the parameters for the WNN generation, see
  [`params_sc_wnn()`](https://gregorlueg.github.io/bixverse/reference/params_sc_wnn.md).

- graph_params:

  Named list. How the fused kNN is turned into the neighbour graph, see
  [`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The `SingleCellsMultiModal` object with the WNN graph and per-cell
modality weights added to `other_data[["wnn"]]`.

## References

Hao et al., Cell, 2021

## Examples

``` r
# a fused RNA + ADT graph, both modalities reduced with PCA first
rna <- generate_single_cell_test_data()
adt <- generate_single_cell_test_data_adt()
dir <- tempfile("bixverse_mm")
dir.create(dir)
object <- load_r_data(
  SingleCellsMultiModal(dir_data = dir),
  counts = rna$counts,
  obs = rna$obs,
  var = rna$var,
  sc_qc_param = params_sc_min_quality(min_unique_genes = 5L),
  .verbose = FALSE
)
object <- add_adt_counts_sc(object, adt_counts = adt$counts, method = "clr")
object <- find_hvg_sc(object, hvg_no = 30L, .verbose = FALSE)
object <- calculate_pca_sc(object, no_pcs = 15L, .verbose = FALSE)
object <- calculate_pca_adt_sc(object, no_pcs = 10L)
object <- generate_wnn_graph_sc(object, .verbose = FALSE)
get_connectivities(object, modality = "wnn")
#> SingleCellConnectivity: 1000 cells, 35379 edges
#>   Method: snn

unlink(dir, recursive = TRUE, force = TRUE)
```
