# Load multiple h5ad files into a single `SingleCells`

Takes a pre-scan result from
[`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md)
and loads all files into a single experiment with global gene QC and
sequential cell indexing.

## Usage

``` r
load_multi_h5ad(
  object,
  prescan_result,
  sc_qc_param = params_sc_min_quality(),
  cell_id_col = NULL,
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- prescan_result:

  Output of
  [`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- cell_id_col:

  Optional string. Column name for cell identifiers in obs.

- streaming:

  Boolean. Use streaming for CSR-to-CSC conversion.

- batch_size:

  Integer. Batch size if `streaming = TRUE`.

- .verbose:

  Boolean.

## Value

The class with updated shape and populated DuckDB.
