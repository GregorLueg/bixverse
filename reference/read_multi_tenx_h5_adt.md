# Read in 10x h5 ADT data from multiple files

Multi-file counterpart to
[`read_tenx_h5_adt()`](https://gregorlueg.github.io/bixverse/reference/read_tenx_h5_adt.md).
Reads the same modality from each input, stacks the cells (rows), and
prefixes each barcode with its `exp_id` so the result matches the
cell_id convention used by
[`load_multi_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_multi_tenx_h5.md).
The feature space is either the intersection or union of features across
inputs; missing features in the union case are filled with zero.

## Usage

``` r
read_multi_tenx_h5_adt(
  h5_paths,
  feature_type = "Antibody Capture",
  gene_universe = c("intersection", "union")
)
```

## Arguments

- h5_paths:

  Character vector of file paths to 10x h5 files. If names are provided,
  these will be used as `exp_id`s; otherwise the file basename is used.

- feature_type:

  String. The feature type to return. Defaults to `"Antibody Capture"`.

- gene_universe:

  One of `"intersection"` or `"union"`.

## Value

A dense matrix of cells x features with `exp_id_barcode` rownames and
feature names as colnames.

## Examples

``` r
# stack the ADT layer of two files, barcodes prefixed by exp_id
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
adt <- generate_single_cell_test_data_adt(
  params_sc_synthetic_data_adt(n_cells = 200L)
)
features <- data.table::data.table(
  id = c(data$var$gene_id, colnames(adt$counts)),
  name = c(data$var$ensembl_id, colnames(adt$counts)),
  feature_type = rep(
    c("Gene Expression", "Antibody Capture"),
    c(ncol(data$counts), ncol(adt$counts))
  )
)
counts <- cbind(data$counts, as(adt$counts, "RsparseMatrix"))
files <- c(a = tempfile(fileext = ".h5"), b = tempfile(fileext = ".h5"))
for (f in files) {
  write_tenx_h5_sc(f, counts, data$obs$cell_id, features)
}
adt_counts <- read_multi_tenx_h5_adt(files)
dim(adt_counts)
#> [1] 400  15
head(rownames(adt_counts), 2)
#> [1] "a_cell_001" "a_cell_002"

unlink(files)
```
