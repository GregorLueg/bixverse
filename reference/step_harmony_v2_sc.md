# Pipeline step: Harmony v2 batch correction

Wraps
[`harmony_v2_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_v2_sc.md)
as an `ScStep`.

## Usage

``` r
step_harmony_v2_sc(
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
