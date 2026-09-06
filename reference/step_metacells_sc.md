# Pipeline step: generate meta cells

Wraps the meta cell generators as an `ScStep`. Unlike the other steps
this one changes the class of the object: it takes a `SingleCells` or
`SingleCellsSubset` and returns a
[`MetaCells()`](https://gregorlueg.github.io/bixverse/reference/MetaCells.md).
Steps that follow it need `MetaCells` methods, which
[`validate_pipeline()`](https://gregorlueg.github.io/bixverse/reference/validate_pipeline.md)
checks up front.

Combined with
[`apply_pipeline_per_group()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline_per_group.md)
this gives you per-group pre-processing (HVG, PCA, batch correction
within a patient, kNN) followed by source-pure meta cells, which you
then hand to
[`merge_meta_cells()`](https://gregorlueg.github.io/bixverse/reference/merge_meta_cells.md).

## Usage

``` r
step_metacells_sc(method = c("bootstrapped", "seacells", "supercells"), ...)
```

## Arguments

- method:

  String. One of `c("bootstrapped", "seacells", "supercells")`.

- ...:

  Arguments passed on to the generator, e.g. `sc_meta_cell_params`,
  `target_size` or `.verbose`.

## Value

An `ScStep`.

## Examples

``` r
# per group pre-processing that ends on source-pure meta cells
sc <- demo_single_cells(prepped = FALSE)
pipeline <- step_hvg_sc(hvg_no = 30L, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = 10L, .verbose = FALSE) %>>%
  step_neighbours_sc(.verbose = FALSE) %>>%
  step_metacells_sc(
    "bootstrapped",
    sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 10L),
    .verbose = FALSE
  )

per_group <- apply_pipeline_per_group(pipeline, sc, group_col = "cell_grp")
merge_meta_cells(per_group, .verbose = FALSE)
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
