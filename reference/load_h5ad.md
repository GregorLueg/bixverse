# Load in h5ad to `single_cell_exp`

This function takes an h5ad file and loads the obs and var data into the
DuckDB of the `single_cell_exp` class and the counts into a
Rust-binarised format for rapid access. During the reading in of the
counts, the log CPM transformation will occur automatically.

## Usage

``` r
load_h5ad(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- h5_path:

  File path to the h5ad object you wish to load in.

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
