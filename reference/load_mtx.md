# Load in mtx/plain text files to `single_cell_exp` (nightly!)

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
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- sc_mtx_io_param:

  List. Please generate this one via
  [`params_sc_mtx_io()`](params_sc_mtx_io.md). Needs to contain:

  - path_mtx - String. Path to the .mtx file

  - path_obs - String. Path to the file containing cell/barcode info.

  - path_var - String. String. Path to the file containing gene/variable
    info.

  - cells_as_rows - Boolean. Do cells represent the rows or columns.

- sc_qc_param:

  List. Output of [`params_sc_min_quality()`](params_sc_min_quality.md).
  A list with the following elements:

  - min_unique_genes - Integer. Minimum number of genes to be detected
    in the cell to be included.

  - min_lib_size - Integer. Minimum library size in the cell to be
    included.

  - min_cells - Integer. Minimum number of cells a gene needs to be
    detected to be included.

  - target_size - Float. Target size to normalise to. Defaults to `1e5`.

- streaming:

  Boolean. Shall the data be streamed during the conversion of CSR to
  CSC. Defaults to `TRUE` and should be used for larger data sets.

- batch_size:

  Integer. If `streaming = TRUE`, how many cells to process in one
  batch. Defaults to `1000L`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
