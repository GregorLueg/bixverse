# Hash of the cell state an artefact must agree with

Order sensitive on purpose.
[`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md)
stores the filter unique and ascending, but memory files written before
0.9.0 can hold an unsorted one; their cached rows follow that order and
must hash as stale.

## Usage

``` r
.sc_state_hash(x)
```

## Arguments

- x:

  The object owning the cache.

## Value

String. The hash of the current cell state.
