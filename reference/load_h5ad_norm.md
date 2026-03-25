# Load in h5ad with normalised counts to `SingleCells`

This function takes an h5ad file where only normalised counts are
available in the X slot and loads the obs and var data into the DuckDB
of the `SingleCells` class and the counts into a Rust-binarised format
for rapid access. Raw counts are reconstructed from the normalised
values using the library sizes stored in a specified obs column.

The reconstruction assumes the normalisation was:
`norm = log1p(x / lib_size * target_size)`

## Usage

``` r
load_h5ad_norm(
  object,
  h5_path,
  obs_lib_size_col,
  target_size,
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  cell_id_col = NULL,
  batch_size = 1000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- h5_path:

  File path to the h5ad object you wish to load in.

- obs_lib_size_col:

  String. Name of the obs column containing the total counts per cell or
  spot (e.g. `"nCount_RNA"`). You will have to check the original obs
  data for this.

- target_size:

  Numeric. The target size used in the original normalisation (e.g.
  `1e4`).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).
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

- cell_id_col:

  Optional string. If a specific column in the h5ad obs data represents
  the cell identifiers, you can specify it here.

- batch_size:

  Integer. If `streaming = TRUE`, how many cells to process in one
  batch. Defaults to `1000L`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
