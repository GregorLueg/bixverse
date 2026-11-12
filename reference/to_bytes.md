# Serialise an object to bytes

Serialises the Rust neighbour graph of a `SingleCellConnectivity` into a
raw vector, which survives
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and qs2 where the
pointer does not. The build method rides along as the attribute
`method`. Read it back with
[`from_bytes()`](https://gregorlueg.github.io/bixverse/reference/from_bytes.md).

## Usage

``` r
to_bytes(x, ...)

# S3 method for class 'SingleCellConnectivity'
to_bytes(x, ...)
```

## Arguments

- x:

  A `SingleCellConnectivity`.

- ...:

  Other parameters.

## Value

A raw vector.

## Examples

``` r
conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
bytes <- to_bytes(conn)
from_bytes(bytes)
#> SingleCellConnectivity: 34 cells, 78 edges
#>   Method: NA
```
