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

- prescan_result:

  Output of
  [`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- cell_id_col:

  Optional string. Column name for cell identifiers in obs.

- csc_mem_gb:

  Optional numeric. Memory in GB for the buffers of the cell-to-gene
  (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL` (default)
  converts in one pass and holds the whole matrix. Set a cap for large
  data sets; every extra phase re-reads the cell file once.

- streaming, batch_size, max_genes_in_memory, cell_batch_size:

  Replaced by `csc_mem_gb` and ignored. **\[deprecated\]**

- .verbose:

  Boolean.

## Value

The class with updated shape and populated DuckDB.

## Examples

``` r
# two files into one experiment, cells tagged by exp_id
files <- c(a = tempfile(fileext = ".h5ad"), b = tempfile(fileext = ".h5ad"))
for (i in seq_along(files)) {
  data <- generate_single_cell_test_data(
    syn_data_params = params_sc_synthetic_data(
      n_cells = 200L,
      n_genes = 40L
    ),
    seed = i
  )
  write_h5ad_sc(files[i], data$counts, data$obs, data$var, .verbose = FALSE)
}
tasks <- prescan_h5ad_files(h5_paths = files, .verbose = FALSE)
dir_data <- tempfile("sc_multi_h5ad")
dir.create(dir_data, recursive = TRUE)
sc <- load_multi_h5ad(
  object = SingleCells(dir_data = dir_data),
  prescan_result = tasks,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  .verbose = FALSE
)
table(sc[["exp_id"]])
#> exp_id
#>   a   b 
#> 200 200 

unlink(c(files, dir_data), recursive = TRUE, force = TRUE)
```
