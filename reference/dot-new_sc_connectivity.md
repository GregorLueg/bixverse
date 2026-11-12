# Wrap a Rust neighbour graph pointer

Low level constructor around the external pointer returned by the Rust
graph builders. Reads the shape back from Rust, so a graph built here
and one rehydrated from bytes end up identical.

## Usage

``` r
.new_sc_connectivity(ptr, method = NA_character_)
```

## Arguments

- ptr:

  External pointer to a Rust `SparseGraph<f32>`.

- method:

  String. How the graph was built, e.g. `"snn"`, `"umap"`, `"knn"`,
  `"bbknn"`. `NA` if unknown.

## Value

A `SingleCellConnectivity` list with

- ptr - The external pointer.

- n - Integer. Number of nodes (cells).

- n_edges - Numeric. Number of undirected edges.

- method - String. How the graph was built.
