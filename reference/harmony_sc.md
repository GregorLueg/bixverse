# Run Harmony

A version of Harmony by Korsunsky et al., implemented in Rust. Performs
batch correction on PCA embeddings and stores the result as a
`"harmony"` embedding in the object.

## Usage

``` r
harmony_sc(
  object,
  batch_column,
  additional_batch_columns = NULL,
  harmony_params = params_sc_harmony(),
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

- harmony_params:

  List. Output of
  [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with a `"harmony"` embedding added. If no PCA embeddings are
found, returns the object unchanged with a warning.
