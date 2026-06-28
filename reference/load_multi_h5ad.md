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
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
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

  Integer. `0L` -\> in-memory, `1L` -\> light streaming, `2L` -\> heavy
  streaming with memory upper boundaries. Defaults to `1L`.

- batch_size:

  Integer. Cell batch size when `streaming = 1L`. Defaults to `1000L`.

- max_genes_in_memory:

  Integer. Maximum genes held in memory at once when `streaming = 2L`.
  Defaults to `2000L`.

- cell_batch_size:

  Integer. Cell batch size when `streaming = 2L`. Defaults to `100000L`.

- .verbose:

  Boolean.

## Value

The class with updated shape and populated DuckDB.
