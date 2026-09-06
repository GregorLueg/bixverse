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
  streaming = 0L,
  .verbose = FALSE
)
table(sc[["exp_id"]])
#> exp_id
#>   a   b 
#> 200 200 

unlink(c(files, dir_data), recursive = TRUE, force = TRUE)
```
