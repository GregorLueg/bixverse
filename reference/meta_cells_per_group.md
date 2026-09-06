# Generate source-pure meta cells and merge them

Splits `object` by `group_col`, optionally runs `pipeline` on each
subset, generates meta cells per group and merges the results into a
single
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
object. This gives you meta cells that never mix cells from two
patients/samples, while still being one object you can run SCENIC,
AUCell or NMF over.

The meta cell generators need an embedding, so `pipeline` will normally
be `step_hvg_sc() %>>% step_pca_sc() %>>% step_neighbours_sc()` unless
every subset already carries one.

## Usage

``` r
meta_cells_per_group(
  object,
  group_col,
  method = c("bootstrapped", "seacells", "supercells"),
  mc_params = list(),
  pipeline = NULL,
  groups = NULL,
  feature_space = c("intersect", "union"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`.

- group_col:

  String. Column in obs used to split.

- method:

  String. One of `c("bootstrapped", "seacells", "supercells")`. Picks
  the meta cell generator.

- mc_params:

  Named list. Arguments passed on to the generator, e.g.
  `list(sc_meta_cell_params = params_sc_bt_metacells(), target_size = 1e5)`.

- pipeline:

  Optional `ScPipeline` applied to each subset before the meta cells are
  generated.

- groups:

  Optional character vector. Restrict to these group values; if `NULL`,
  all unique values of `group_col` are used.

- feature_space:

  String. One of `c("intersect", "union")`. Passed to
  [`merge_meta_cells()`](https://gregorlueg.github.io/bixverse/reference/merge_meta_cells.md).
  Irrelevant here as all groups share the gene space of the parent
  object.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A merged
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
object with a `source_id` column in its observation table.

## Examples

``` r
# meta cells that never mix two cell groups
sc <- demo_single_cells(prepped = FALSE)
prep <- step_hvg_sc(hvg_no = 30L, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = 10L, .verbose = FALSE) %>>%
  step_neighbours_sc(.verbose = FALSE)
meta_cells_per_group(
  object = sc,
  group_col = "cell_grp",
  method = "bootstrapped",
  mc_params = list(
    sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 10L),
    .verbose = FALSE
  ),
  pipeline = prep,
  .verbose = FALSE
)
#> Single cell experiment (Meta Cells).
#>   Meta cell method: meta_cells_hdwgcna
#>   Merged: TRUE
#>   No meta cells: 30
#>   No genes: 50
#>   No cells aggregated: 297
#>   No obs rows in source: 500
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
