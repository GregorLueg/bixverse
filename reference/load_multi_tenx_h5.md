# Load multiple 10x CellRanger h5 files into a single `SingleCells`

Takes the result of
[`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md)
and loads all inputs into a single experiment with global gene QC and
sequential cell indexing. The feature space is determined by the prescan
(intersection or union of gene ids).

## Usage

``` r
load_multi_tenx_h5(
  object,
  prescan_result,
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

- prescan_result:

  Output of
  [`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md).

- sc_qc_param:

  List. Output of
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

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
# two 10x h5 files into one experiment
files <- c(a = tempfile(fileext = ".h5"), b = tempfile(fileext = ".h5"))
for (i in seq_along(files)) {
  data <- generate_single_cell_test_data(
    syn_data_params = params_sc_synthetic_data(
      n_cells = 200L,
      n_genes = 40L
    ),
    seed = i
  )
  write_tenx_h5_sc(
    f_path = files[i],
    counts = data$counts,
    barcodes = data$obs$cell_id,
    features = data.table::data.table(
      id = data$var$gene_id,
      name = data$var$ensembl_id,
      feature_type = "Gene Expression"
    )
  )
}
scan_res <- prescan_tenx_h5_files(h5_paths = files, .verbose = FALSE)

dir_data <- tempfile("sc_multi_tenx")
dir.create(dir_data, recursive = TRUE)
sc <- load_multi_tenx_h5(
  object = SingleCells(dir_data = dir_data),
  prescan_result = scan_res,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 5L,
    min_lib_size = 25L,
    min_cells = 5L
  ),
  .verbose = FALSE
)
dim(sc)
#> [1] 400  40

unlink(c(files, dir_data), recursive = TRUE, force = TRUE)
```
