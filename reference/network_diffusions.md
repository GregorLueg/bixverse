# Network diffusion class (deprecated)

**\[deprecated\]**

This constructor has been renamed to
[`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

## Usage

``` r
network_diffusions(edge_data_frame, weighted, directed)
```

## Arguments

- edge_data_frame:

  data.table that contains the edge information. It is expected to have
  the columns 'from' and 'to'.

- weighted:

  Boolean. Is the graph weighted. If set to TRUE, the `edge_data_frame`
  needs to have a weight column.

- directed:

  Boolean. Shall the graph be stored as directed.

## Value

Returns a
[`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md)
object.
