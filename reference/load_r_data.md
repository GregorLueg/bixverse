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
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
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

- csc_mem_gb:

  Optional numeric. Memory in GB for the buffers of the cell-to-gene
  (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL` (default)
  converts in one pass and holds the whole matrix. Set a cap for large
  data sets; every extra phase re-reads the cell file once.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Replaced by `csc_mem_gb` and ignored. **\[deprecated\]**

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

It will populate the files on disk and return the class with updated
shape information.

## Examples

``` r
# straight from a dgRMatrix in memory onto disk
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
dir_data <- tempfile("sc_r_data")
dir.create(dir_data, recursive = TRUE)
sc <- load_r_data(
  object = SingleCells(dir_data = dir_data),
  counts = data$counts,
  obs = data$obs,
  var = data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  .verbose = FALSE
)
dim(sc)
#> [1] 200  40

unlink(dir_data, recursive = TRUE, force = TRUE)
```
