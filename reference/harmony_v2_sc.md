# Run Harmony v2

A version of Harmony v2 by Patikas et al., 2026, implemented in Rust.
Performs batch correction on PCA embeddings and stores the result as a
`"harmony_v2"` embedding in the object.

## Usage

``` r
harmony_v2_sc(
  object,
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony_v2(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

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

The object with a `"harmony_v2"` embedding added. If no PCA embeddings
are found, returns the object unchanged with a warning.
