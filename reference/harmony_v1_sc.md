# Run Harmony v1

The original Harmony by Korsunsky et al., implemented in Rust. Performs
batch correction on PCA embeddings and stores the result as a
`"harmony_v1"` embedding in the object.
[`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md)
runs Harmony v2, which is the default.

## Usage

``` r
harmony_v1_sc(
  object,
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

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
  [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md).

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The object with a `"harmony_v1"` embedding added. If no PCA embeddings
are found, returns the object unchanged with a warning.

## Examples

``` r
# the original Harmony, writing its own embedding
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
sc <- harmony_v1_sc(sc, batch_column = "batch_index", .verbose = FALSE)
dim(get_embedding(sc, "harmony_v1"))
#> [1] 600  10

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
