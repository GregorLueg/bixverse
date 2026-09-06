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
  streaming = 1L,
  cell_id_col = NULL,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
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
  spot (e.g. `"nCount_RNA"`).

- target_size:

  Numeric. The target size used in the original normalisation (e.g.
  `1e4`).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- streaming:

  Integer. `0L` -\> in-memory, `1L` -\> light streaming, `2L` -\> heavy
  streaming with memory upper boundaries. Controls memory pressure
  during CSR-to-CSC conversion. Defaults to `1L`.

- cell_id_col:

  Optional string. If a specific column in the h5ad obs data represents
  the cell identifiers, you can specify it here.

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

It will populate the files on disk and return the class with updated
shape information.

## Examples

``` r
# an h5ad holding log1p(x / lib_size * 1e4), raw counts reconstructed on read
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
lib_size <- Matrix::rowSums(data$counts)
norm_counts <- data$counts
norm_counts@x <- log1p(
  norm_counts@x / rep(lib_size, diff(norm_counts@p)) * 1e4
)
obs <- data.table::copy(data$obs)[, total_counts := lib_size]

f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, norm_counts, obs, data$var, .verbose = FALSE)
dir_data <- tempfile("sc_h5ad_norm")
dir.create(dir_data, recursive = TRUE)
sc <- load_h5ad_norm(
  object = SingleCells(dir_data = dir_data),
  h5_path = f_path,
  obs_lib_size_col = "total_counts",
  target_size = 1e4,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  streaming = 0L,
  .verbose = FALSE
)
dim(sc)
#> [1] 200  40

unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
```
