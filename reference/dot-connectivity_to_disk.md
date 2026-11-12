# Swap the Rust pointer of a cached graph for bytes

Used by
[`save_sc_exp_to_disk()`](https://gregorlueg.github.io/bixverse/reference/save_sc_exp_to_disk.md).
Rewrites the elements of the list only, so the provenance stamp
attribute is kept. A graph already in its on-disk form is returned as
is.

## Usage

``` r
.connectivity_to_disk(x)
```

## Arguments

- x:

  The cached graph, or `NULL`.

## Value

The graph with `bytes` in place of `ptr`.
