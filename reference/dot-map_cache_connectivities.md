# Apply a function to the neighbour graph of a cache, if there is one

Used to swap the Rust graph for bytes on save and back on load. An
absent graph is left absent rather than written as `NULL`, which on a
list would delete the slot.

## Usage

``` r
.map_cache_connectivities(x, f)
```

## Arguments

- x:

  The cache list (an `ScCache` or the `wnn` slot), or `NULL`.

- f:

  Function. Takes and returns the graph.

## Value

The cache with `f` applied to `connectivities`.
