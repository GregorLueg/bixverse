# Stream in h5ad to `SingleCells` (alias)

Convenience alias for `load_h5ad(streaming = 2L)`. Kept for backwards
compatibility - forwards directly to
[`load_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad.md)
with heavy streaming enabled. Prefer calling `load_h5ad` directly with
an explicit `streaming` level.

## Usage

``` r
stream_h5ad(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- h5_path:

  File path to the h5ad object.

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- raw_count_slot:

  Where raw counts live. `"auto"` detects per file via
  [`detect_raw_count_slot()`](https://gregorlueg.github.io/bixverse/reference/detect_raw_count_slot.md);
  otherwise one of `"X"`, `"raw.X"`, `"layers.counts"`.

- max_genes_in_memory:

  Integer. Genes held in memory at once. Defaults to `2000L`.

- cell_batch_size:

  Integer. Cell batch size. Defaults to `100000L`.

- .verbose:

  Boolean.

## Value

The class with updated shape information.

## Examples

``` r
# same as load_h5ad(streaming = 2L)
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
dir_data <- tempfile("sc_stream")
dir.create(dir_data, recursive = TRUE)
sc <- stream_h5ad(
  object = SingleCells(dir_data = dir_data),
  h5_path = f_path,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  .verbose = FALSE
)
dim(sc)
#> [1] 200  40

unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
```
