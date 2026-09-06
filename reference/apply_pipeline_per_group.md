# Apply a pipeline independently to each group of a `SingleCells` object

Splits `object` by `group_col`, applies `pipeline` to each subset, and
returns a named list of processed `SingleCellsSubset`s. Useful for
per-sample / per-cell-type re-analysis where the same chain (HVG, PCA,
neighbours, clusters, ...) is run on each group, e.g. sample-pure
metacell generation followed by an external merge.

## Usage

``` r
apply_pipeline_per_group(
  pipeline,
  object,
  group_col,
  groups = NULL,
  progress = FALSE
)
```

## Arguments

- pipeline:

  `ScPipeline`.

- object:

  `SingleCells`.

- group_col:

  String. Column in obs used to split.

- groups:

  Optional character vector. Restrict to these group values; if `NULL`,
  all unique values of `group_col` are used.

- progress:

  Boolean. Shall big progress messages be printed to the console.
  Defaults to `FALSE`.

## Value

Named list of processed objects, names being the group values. Usually
`SingleCellsSubset`, or `MetaCells` if the pipeline ends on
[`step_metacells_sc()`](https://gregorlueg.github.io/bixverse/reference/step_metacells_sc.md),
in which case
[`merge_meta_cells()`](https://gregorlueg.github.io/bixverse/reference/merge_meta_cells.md)
puts them back together.

## Examples

``` r
# the same chain re-run inside each cell type
sc <- demo_single_cells(prepped = FALSE)
p <- sc_pipeline() %>>% step_hvg_sc(hvg_no = 20L, .verbose = FALSE)
res <- apply_pipeline_per_group(p, sc, group_col = "cell_grp")
names(res)
#> [1] "cell_type_1" "cell_type_2" "cell_type_3"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
