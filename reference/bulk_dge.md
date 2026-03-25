# Bulk RNAseq differential gene expression class (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

## Usage

``` r
bulk_dge(
  raw_counts,
  meta_data,
  variable_info = NULL,
  alternative_gene_id = NULL
)
```

## Arguments

- raw_counts:

  matrix. The raw count matrix. Rows = genes, columns = samples.

- meta_data:

  data.table. Metadata information on the samples. Expects to have a
  `sample_id` column.

- variable_info:

  data.table. Metadata information on the features. Defaults to `NULL`.

- alternative_gene_id:

  String. Optional alternative gene identifier to be used. Must be a
  column of `variable_info`.

## Value

Returns a
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
object.
