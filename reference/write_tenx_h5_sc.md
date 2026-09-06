# Helper function to write data to a 10x CellRanger-style h5 file

Helper function to write data to a 10x CellRanger-style h5 file

## Usage

``` r
write_tenx_h5_sc(
  f_path,
  counts,
  barcodes,
  features,
  version = c("v3", "v2"),
  overwrite = TRUE
)
```

## Arguments

- f_path:

  String. Output path.

- counts:

  Sparse matrix (`dgRMatrix` or `dgCMatrix`), cells x features.

- barcodes:

  Character. Cell barcodes, length `nrow(counts)`.

- features:

  data.table with `id` and `name` of length `ncol(counts)`. For v3 may
  include `feature_type`; defaults to `"Gene Expression"` if absent.

- version:

  One of `"v3"` or `"v2"`.

- overwrite:

  Boolean.

## Value

Invisible.

## Examples

``` r
# a CellRanger v3 style h5, ready for load_tenx_h5()
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5")
write_tenx_h5_sc(
  f_path = f_path,
  counts = data$counts,
  barcodes = data$obs$cell_id,
  features = data.table::data.table(
    id = data$var$gene_id,
    name = data$var$ensembl_id
  )
)
read_tenx_h5_metadata(f_path)$dims
#> obs var 
#> 200  40 

unlink(f_path)
```
