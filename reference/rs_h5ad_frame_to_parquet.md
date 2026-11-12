# Write the obs or var dataframe of an h5ad file to Parquet

**\[experimental\]** Writes the table the obs / var DuckDB tables are
built from: `cell_idx` / `gene_idx` (1-based after filtering), `cell_id`
/ `gene_id`, then the remaining columns snake_cased. Categoricals are
written as strings; their levels come back for the caller to cast to an
ENUM.

## Usage

``` r
rs_h5ad_frame_to_parquet(h5_path, frame, keep, id_col, out_path)
```

## Arguments

- h5_path:

  String. Path to the h5ad file.

- frame:

  String. One of `c("obs", "var")`.

- keep:

  Optional integer vector. 0-indexed(!) rows to keep, in output order.
  `NULL` keeps all.

- id_col:

  Optional string. Column to use as identifier when the group has no
  index; the first column otherwise.

- out_path:

  String. Parquet file to write.

## Value

A list with

- names - Character vector. Column names in file order.

- types - Character vector. One of `"integer"`, `"double"`, `"boolean"`,
  `"string"` or `"categorical"` per column.

- categories - List. The levels in order for a categorical column,
  `NULL` otherwise.

- skipped - Character vector. Source columns left out for an unsupported
  encoding.

- n_rows - Numeric. Number of rows written.
