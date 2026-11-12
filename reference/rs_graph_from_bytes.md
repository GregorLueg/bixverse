# Read a Rust neighbour graph back from bytes

**\[experimental\]** Inverse of
[`rs_graph_to_bytes()`](https://gregorlueg.github.io/bixverse/reference/rs_graph_to_bytes.md).
The CSR is validated on the way in.

## Usage

``` r
rs_graph_from_bytes(bytes)
```

## Arguments

- bytes:

  Raw vector written by
  [`rs_graph_to_bytes()`](https://gregorlueg.github.io/bixverse/reference/rs_graph_to_bytes.md).

## Value

External pointer to the Rust graph.
