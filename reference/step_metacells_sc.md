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
if (FALSE) { # \dontrun{
pipeline <- step_hvg_sc() %>>%
  step_pca_sc(no_pcs = 20L) %>>%
  step_harmony_sc(batch_column = "plate") %>>%
  step_neighbours_sc(embd_to_use = "harmony") %>>%
  step_metacells_sc("bootstrapped")

per_patient <- apply_pipeline_per_group(pipeline, sc_obj, "patient_id")
mc <- merge_meta_cells(per_patient)
} # }
```
