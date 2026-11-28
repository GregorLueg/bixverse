# Open an experiment's gene activity store

Opens the store
[`gene_activity_sc()`](https://gregorlueg.github.io/bixverse/reference/gene_activity_sc.md)
wrote, with the cell filter of `object`.

## Usage

``` r
get_gene_activity(object, dir_data = NULL)
```

## Arguments

- object:

  `SingleCells` class the gene activity was counted for.

- dir_data:

  Optional string. Where it was written; `NULL` for the default
  `<dir_data>/gene_activity`.

## Value

The gene activity experiment as a `SingleCells`. Errors if there is
none, or if it was counted for other cells (a re-ingest since).
