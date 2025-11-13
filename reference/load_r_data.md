# Load in data directly from R objects.

This function loads in data directly from R objects.

## Usage

``` r
load_r_data(
  object,
  counts,
  obs,
  var,
  sc_qc_param = params_sc_min_quality(),
  batch_size = 1000L,
  streaming = TRUE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- counts:

  Sparse matrix. The cells represent the rows, the genes the indices.
  Needs to be `"dgRMatrix"`.

- obs:

  data.table. The data.table representing the observations, i.e., cell
  information.

- var:

  data.table. The data.table representing the features, i.e., the
  feature information.

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

- batch_size:

  Integer. If `streaming = TRUE`, how many cells to process in one
  batch. Defaults to `1000L`.

- streaming:

  Boolean. Shall the data be streamed during the conversion of CSR to
  CSC. Defaults to `TRUE` and should be used for larger data sets.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
