# Calculate VISION scores (with auto-correlation scores)

Calculates VISION-type scores for pathways based on DeTomaso, et al.
Compared to other score types, you can also calculate delta-type scores
between positive and negative gene indices, think epithelial vs
mesenchymal gene signature, etc. Additionally, this function also
calculates the auto- correlation values, answering the question if a
given signature shows non- random enrichment on the kNN graph. The kNN
graph (and distance measures) will be generated on-the-fly based on the
embedding you wish to use.

## Usage

``` r
vision_w_autocor_sc(
  object,
  gs_list,
  embd_to_use,
  no_embd_to_use = NULL,
  use_knn = TRUE,
  vision_params = params_sc_vision(),
  streaming = NULL,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- gs_list:

  Named nested list. Every element must itself be a list with at least a
  `"pos"` element holding the gene identifiers, and optionally a `"neg"`
  one. A bare character vector is not accepted. The gene identifiers
  need to be part of the variables of the object.

- embd_to_use:

  String. The embedding to use. Whichever you chose, it needs to be part
  of the object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- use_knn:

  Boolean. Shall the internal kNN be used. If set to yes, you need to
  ensure consistency.

- vision_params:

  List with vision parameters, see
  [`params_sc_vision()`](https://gregorlueg.github.io/bixverse/reference/params_sc_vision.md)
  with the following elements:

  - n_perm - Integer. Number of random permutations

  - n_cluster - Integer. Number of random clusters to generate to
    associate each set with.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect. Ignored when applied to
  `MetaCells`.

- random_seed:

  Integer. The random seed.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A list with the following elements:

- vision_matrix - Matrix of cells x signatures with the VISION pathway
  scores as values.

- auto_cor_dt - data.table with the auto-correlation results per gene
  set, i.e., `auto_cor` (1 - Gaery's C), `p_val` and `fdr`.

## References

DeTomaso, et al., Nat. Commun., 2019

## Examples

``` r
# scores plus whether they sit non-randomly on the kNN graph
sc <- demo_single_cells()
gs_list <- list(
  programme_a = list(pos = get_gene_names(sc)[1:10]),
  programme_b = list(pos = get_gene_names(sc)[21:30])
)
res <- vision_w_autocor_sc(
  sc,
  gs_list = gs_list,
  embd_to_use = "pca",
  vision_params = params_sc_vision(n_perm = 50L, n_cluster = 3L),
  .verbose = FALSE
)
res$auto_cor_dt
#>    gene_set_name  auto_cor      p_val        fdr
#>           <char>     <num>      <num>      <num>
#> 1:   programme_a 0.7248461 0.01960784 0.01960784
#> 2:   programme_b 0.7437420 0.01960784 0.01960784

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
