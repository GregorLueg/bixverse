# Load multiple 10x CellRanger h5 files into a single `SingleCells`

Takes the result of
[`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md)
and loads all inputs into a single experiment with global gene QC and
sequential cell indexing. The feature space is determined by the prescan
(intersection or union of gene ids).

## Usage

``` r
load_multi_tenx_h5(
  object,
  prescan_result,
  sc_qc_param = params_sc_min_quality(),
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
  [`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- streaming:

  Integer. CSR-to-CSC conversion mode. `0L` -\> in-memory, `1L` -\>
  light streaming, `2L` -\> heavy streaming with memory upper
  boundaries. Defaults to `1L`.

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
