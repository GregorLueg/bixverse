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
