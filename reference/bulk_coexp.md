# Bulk RNAseq co-expression modules (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Usage

``` r
bulk_coexp(raw_data, meta_data, variable_info = NULL)
```

## Arguments

- raw_data:

  The raw count matrix. Rows = samples, columns = features.

- meta_data:

  data.table. Metadata information on the samples. Expects to have a
  `sample_id` column.

- variable_info:

  data.table. Metadata information on the features. This is an optional
  table.

## Value

Returns a
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
object.
