# Extract the CSR adjacency of a Rust neighbour graph

**\[experimental\]** The symmetric adjacency, rows sorted, as ScanPy
stores it under `obsp/connectivities`.

## Usage

``` r
rs_graph_to_csr(ptr)
```

## Arguments

- ptr:

  External pointer to the Rust graph.

## Value

A list with

- indptr - Integer vector. Row pointers (0-indexed).

- indices - Integer vector. Column indices (0-indexed).

- data - Numeric vector. Edge weights.
