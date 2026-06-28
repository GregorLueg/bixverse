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
kernelised pseudo-distance metric), the sNN graph as an igraph and a
table of per-cell modality weights. This graph can subsequently be used
for clustering and 2D embedding via the `"wnn"` modality.

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
  full_snn = TRUE,
  pruning = 1/15,
  snn_similarity = "jaccard",
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

- full_snn:

  Boolean. Shall the full shared nearest neighbour graph be generated
  that generates edges between all cells instead of between only
  neighbours.

- pruning:

  Numeric. Weights below this threshold will be set to 0 in the
  generation of the sNN graph. Seurat uses for example 1/15 with k = 20.

- snn_similarity:

  String. One of `c("rank", "jaccard")`. The Jaccard similarity
  calculates the Jaccard between the neighbours, whereas the rank method
  calculates edge weights based on the ranking of shared neighbours. For
  the rank method, the weight is determined by finding the shared
  neighbour with the lowest combined rank across both cells, where
  lower-ranked (closer) shared neighbours result in higher edge weights
  Both methods produce weights normalised to the range `⁠[0, 1]`⁠.

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
