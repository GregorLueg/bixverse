# Stamp the cell IDs of a freshly written store into its headers

Only their hash is stored. Stores that share a cell axis, such as a gene
activity store and its source, carry the same one, which is how
[`get_gene_activity()`](https://gregorlueg.github.io/bixverse/reference/get_gene_activity.md)
detects a re-ingest.

## Usage

``` r
.stamp_store_cell_ids(dir_data)
```

## Arguments

- dir_data:

  String. The directory with the store and its DuckDB.

## Value

Invisible `NULL`.
