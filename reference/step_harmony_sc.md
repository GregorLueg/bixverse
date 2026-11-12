# Pipeline step: Harmony batch correction

Wraps
[`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md)
(Harmony v2) as an `ScStep`.

## Usage

``` r
step_harmony_sc(
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony_v2(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- batch_column:

  String. Column name in the object containing the primary batch labels.

- additional_batch_columns:

  Optional character vector. Additional batch columns to regress out. If
  `NULL`, only the primary batch column is used.

- modality:

  String. One of `c("rna", "adt")`. You can only use `"adt"` on
  `SingleCellsMultiModal` class.

- harmony_params:

  List. Output of
  [`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

An `ScStep`.

## Examples

``` r
# correct the PCA, then build the graph on the corrected embedding
step_harmony_sc(batch_column = "batch_index") %>>%
  step_neighbours_sc(embd_to_use = "harmony")
#> <ScPipeline> 2 steps
#>   1. harmony     batch_column = "batch_index", additional_batch_columns = NULL, modality = c("rna", "adt"), harmony_params = <params_sc_harmony_v2>, seed = 42L, .verbose = TRUE
#>   2. neighbours  embd_to_use = "harmony", no_embd_to_use = NULL, modality = c("rna", "adt"), neighbours_params = <list>, seed = 42L, .verbose = TRUE
```
