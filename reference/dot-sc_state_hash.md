# Hash of the cell state an artefact must agree with

Order sensitive on purpose.
[`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md)
does not sort, and the PCA rows come back from Rust in exactly
`cell_indices` order, so row order is what aligns a cached matrix to a
cell list.

## Usage

``` r
.sc_state_hash(x)
```

## Arguments

- x:

  The object owning the cache.

## Value

String. The hash of the current cell state.
