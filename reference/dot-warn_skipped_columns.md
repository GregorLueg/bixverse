# Warn about the metadata columns the Parquet writer skipped

Warn about the metadata columns the Parquet writer skipped

## Usage

``` r
.warn_skipped_columns(info, group)
```

## Arguments

- info:

  List. Output of
  [`rs_h5ad_frame_to_parquet()`](https://gregorlueg.github.io/bixverse/reference/rs_h5ad_frame_to_parquet.md).

- group:

  String. The h5ad group, e.g. `"/obs"`.

## Value

Invisible `NULL`. Warns once per skipped column.
