# Serialise a Rust neighbour graph to bytes

**\[experimental\]** The byte layout carries a magic, a format version
and the float width, so reading back a foreign or outdated blob errors
rather than decoding into garbage.

## Usage

``` r
rs_graph_to_bytes(ptr)
```

## Arguments

- ptr:

  External pointer to the Rust graph.

## Value

Raw vector.
