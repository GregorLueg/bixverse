# Merge meta cell objects into one

Row-binds several
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md)
objects into a single one. The typical use case is sample-pure meta
cells: generate meta cells per patient (see
[`meta_cells_per_group()`](https://gregorlueg.github.io/bixverse/reference/meta_cells_per_group.md)),
then merge them so that methods like SCENIC, AUCell or NMF can run
across the full set.

The normalised counts are carried over as generated. They are normalised
per meta cell, so row-binding leaves them valid. Caches (PCA, kNN, sNN,
embeddings) are per source and are dropped; recompute them on the merged
object.

## Usage

``` r
merge_meta_cells(
  inputs,
  source_ids = NULL,
  feature_space = c("intersect", "union"),
  prefix_ids = TRUE,
  .verbose = TRUE
)
```

## Arguments

- inputs:

  List of `MetaCells` objects.

- source_ids:

  Optional character vector of the same length as `inputs` with the
  source (e.g. patient) identifiers. Defaults to `names(inputs)` and
  falls back to `source_01`, `source_02`, ... Needs to be unique.

- feature_space:

  String. One of `c("intersect", "union")`. Controls how differing gene
  spaces are resolved. With `"union"` genes missing from an input become
  structural zeros for its meta cells. Irrelevant when all inputs came
  from the same source object, as their gene spaces are then identical.

- prefix_ids:

  Boolean. Prefix the meta cell identifiers with the source identifier.
  If `FALSE`, duplicated identifiers across inputs are an error.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A single `MetaCells` object with all meta cells of the inputs and
`is_merged` set to `TRUE`. The observation table gains a `source_id`
column; `original_cell_idx` stays in the index space of its own source,
which is why methods that resolve it against the source single cell data
([`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md),
[`calc_manifold_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_manifold_metrics.md))
refuse to run on the result. `other_data` holds the source identifiers
under `sources`.

## Examples

``` r
# \donttest{
# meta cells generated per source, then pooled for downstream methods
sc_1 <- demo_single_cells(seed = 1L)
sc_2 <- demo_single_cells(seed = 2L)
params <- params_sc_bt_metacells(target_no_metacells = 25L)
merged <- merge_meta_cells(
  list(
    donor_1 = generate_bt_meta_cells_sc(sc_1, params, .verbose = FALSE),
    donor_2 = generate_bt_meta_cells_sc(sc_2, params, .verbose = FALSE)
  ),
  .verbose = FALSE
)
merged
#> Single cell experiment (Meta Cells).
#>   Meta cell method: meta_cells_hdwgcna
#>   Merged: TRUE
#>   No meta cells: 50
#>   No genes: 50
#>   No cells aggregated: 379
#>   No obs rows in source: 500
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
#>   Stale artefacts: none

unlink(c(sc_1@dir_data, sc_2@dir_data), recursive = TRUE, force = TRUE)
# }
```
