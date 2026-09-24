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
get_snn_graph(object, modality = "wnn")
#> IGRAPH 36076bb U-W- 1000 35379 -- 
#> + attr: weight (e/n)
#> + edges from 36076bb:
#>  [1]  1-- 4  5-- 8  3-- 9  1--13  4--13 13--16  1--16 14--17  3--18 10--19
#> [11] 14--20 17--20  5--20 12--21  1--22 16--22  5--23 20--23  8--23 15--24
#> [21] 12--24  9--24 21--24  4--25  5--26 20--26 23--26 17--26 12--27 21--27
#> [31] 24--27  6--27 25--28  2--29 12--30 15--30 24--30 28--31 17--32  9--33
#> [41] 15--33 31--34  4--34 32--35 11--35 15--36 27--36 30--36 24--36  1--37
#> [51] 13--37 16--37 22--37 28--37 14--38 17--38 20--38 27--39 33--39  9--39
#> [61] 24--39 12--39 36--39 15--39 28--40 31--40  2--41 29--41 11--41 35--41
#> [71] 18--42  3--42  9--42 22--43 10--43  1--43 19--43 32--44 35--44 17--44
#> + ... omitted several edges

unlink(dir, recursive = TRUE, force = TRUE)
```
