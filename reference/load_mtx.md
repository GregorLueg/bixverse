# Load in mtx/plain text files to `SingleCells`

This is a helper function to load in mtx files and corresponding plain
text files. It will automatically filter out low quality cells and only
keep high quality cells. Under the hood DucKDB and high performance Rust
binary files are being used to store the counts.

## Usage

``` r
load_mtx(
  object,
  sc_mtx_io_param = params_sc_mtx_io(),
  sc_qc_param = params_sc_min_quality(),
  mtx_streaming = TRUE,
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

- sc_mtx_io_param:

  List. Please generate this one via
  [`params_sc_mtx_io()`](https://gregorlueg.github.io/bixverse/reference/params_sc_mtx_io.md).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- mtx_streaming:

  Boolean. Shall the .mtx file ingestion itself be streamed (via
  temp-file bucketing). Recommended for large mtx files. Defaults to
  `TRUE`.

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

The class with updated shape information.
