# Write the barcodes or features of a 10x h5 file to Parquet

**\[experimental\]** Obs is `cell_idx`, `cell_id`. Var is `gene_idx`,
`gene_id`, `gene_name`, plus `feature_type` for v3. Features are read in
file order across all modalities, so `keep` indexes into that order.

## Usage

``` r
rs_tenx_frame_to_parquet(h5_path, version, frame, keep, out_path)
```

## Arguments

- h5_path:

  String. Path to the 10x h5 file.

- version:

  String. One of `c("v2", "v3")`.

- frame:

  String. One of `c("obs", "var")`.

- keep:

  Optional integer vector. 0-indexed(!) rows to keep, in output order.
  `NULL` keeps all.

- out_path:

  String. Parquet file to write.

## Value

See
[`rs_h5ad_frame_to_parquet()`](https://gregorlueg.github.io/bixverse/reference/rs_h5ad_frame_to_parquet.md).
