# Load in h5ad to `SingleCells`

This function takes an h5ad file and loads the obs and var data into the
DuckDB of the `SingleCells` class and the counts into a Rust-binarised
format for rapid access. During the reading in of the counts, the log
CPM transformation will occur automatically.

## Usage

``` r
load_h5ad(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  cell_id_col = NULL,
  h5ad_streaming = TRUE,
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

- h5_path:

  File path to the h5ad object you wish to load in.

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

- raw_count_slot:

  Where raw counts live. `"auto"` detects per file via
  [`detect_raw_count_slot()`](https://gregorlueg.github.io/bixverse/reference/detect_raw_count_slot.md);
  otherwise one of `"X"`, `"raw.X"`, `"layers.counts"`.

- cell_id_col:

  Optional string. If a specific column in the h5ad obs data is
  representing the cell identifiers, you can specify it here.

- h5ad_streaming:

  Boolean. Stream the h5ad counts into the cell-based binary in batches
  instead of materialising the filtered matrix first. Recommended for
  large files. Defaults to `TRUE`.

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
# round trip through a sparse h5ad file
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
dir_data <- tempfile("sc_h5ad")
dir.create(dir_data, recursive = TRUE)
sc <- load_h5ad(
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
