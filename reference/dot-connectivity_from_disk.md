# Rehydrate a cached graph loaded off disk

Used by
[`load_existing()`](https://gregorlueg.github.io/bixverse/reference/load_existing.md),
through the cache migrations. A graph in its on-disk form gets its
pointer back. A legacy igraph is converted into a
`SingleCellConnectivity` and dropped. The provenance stamp is kept in
both cases.

## Usage

``` r
.connectivity_from_disk(x)
```

## Arguments

- x:

  The cached graph, or `NULL`.

## Value

A `SingleCellConnectivity`, or `x` as is.
