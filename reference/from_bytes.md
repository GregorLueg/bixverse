# Read a SingleCellConnectivity back from bytes

Inverse of
[`to_bytes()`](https://gregorlueg.github.io/bixverse/reference/to_bytes.md).
Bytes from a different structure, an older byte layout or another float
width error rather than decode into garbage.

## Usage

``` r
from_bytes(bytes)
```

## Arguments

- bytes:

  Raw vector written by
  [`to_bytes()`](https://gregorlueg.github.io/bixverse/reference/to_bytes.md).

## Value

A `SingleCellConnectivity`.

## Examples

``` r
conn <- as_sc_connectivity(igraph::make_graph("Zachary"))
from_bytes(to_bytes(conn))
#> SingleCellConnectivity: 34 cells, 78 edges
#>   Method: NA
```
