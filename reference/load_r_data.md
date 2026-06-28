# Load in data directly from R objects.

This function loads in data directly from R objects. The counts matrix
must be a `dgRMatrix` (rows = cells, columns = genes).

## Usage

``` r
load_r_data(
  object,
  counts,
  obs,
  var,
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

- counts:

  Sparse matrix. The cells represent the rows, the genes the columns.
  Needs to be a `"dgRMatrix"`.

- obs:

  data.table. Cell metadata. Must have one row per cell in the same
  order as `counts`.

- var:

  data.table. Feature metadata. Must have one row per gene in the same
  order as `counts`.

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

  Integer. CSR-to-CSC conversion mode. `0L` -\> in-memory (fastest,
  highest memory), `1L` -\> light streaming with cell batching, `2L` -\>
  heavy streaming with memory upper boundaries. Defaults to `1L`.

- batch_size:

  Integer. Cell batch size when `streaming = 1L`. Defaults to `1000L`.

- max_genes_in_memory:

  Integer. Maximum genes held in memory at once when `streaming = 2L`.
  Defaults to `2000L`.

- cell_batch_size:

  Integer. Cell batch size when `streaming = 2L`. Defaults to `100000L`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.
