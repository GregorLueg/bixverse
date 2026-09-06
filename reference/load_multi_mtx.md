# Load multiple mtx directories into a single `SingleCells`

Takes the result of
[`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md)
and loads all inputs into a single experiment with global gene QC and
sequential cell indexing. The feature space is the **intersection** of
input gene IDs.

## Usage

``` r
load_multi_mtx(
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
  [`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md).

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

## Examples

``` r
# two CellRanger directories into one experiment
dirs <- c(tempfile("cr_a"), tempfile("cr_b"))
for (i in seq_along(dirs)) {
  dir.create(dirs[i], recursive = TRUE)
  data <- generate_single_cell_test_data(
    syn_data_params = params_sc_synthetic_data(
      n_cells = 200L,
      n_genes = 40L
    ),
    seed = i
  )
  write_cellranger_output(
    dirs[i], data$counts, data$obs, data$var,
    rows = "cells", format_type = "csv", .verbose = FALSE
  )
}
scan_res <- prescan_mtx_dirs(
  dirs = dirs,
  exp_ids = c("a", "b"),
  cells_as_rows = TRUE,
  has_hdr = TRUE,
  .verbose = FALSE
)

dir_data <- tempfile("sc_multi_mtx")
dir.create(dir_data, recursive = TRUE)
sc <- load_multi_mtx(
  object = SingleCells(dir_data = dir_data),
  prescan_result = scan_res,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  streaming = 0L,
  .verbose = FALSE
)
dim(sc)
#> [1] 400  40

unlink(
  c(dirs, dir_data, scan_res$temp_files),
  recursive = TRUE,
  force = TRUE
)
```
