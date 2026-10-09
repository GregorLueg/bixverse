# Recompute the layout of a Bonsai tree

Lays the finished tree out again without searching. The tree is
renumbered on the Rust side, so the inferred ancestors can come back
with different node numbers; the leaves, and so the cells, keep theirs.

## Usage

``` r
relayout_bonsai(
  x,
  layout = c("equal_angle", "equal_daylight", "dendrogram"),
  hyperbolic = FALSE
)
```

## Arguments

- x:

  A `BonsaiTree` object.

- layout:

  String. One of `c("equal_angle", "equal_daylight", "dendrogram")`.

- hyperbolic:

  Boolean. Project the layout onto the hyperbolic disk.

## Value

The `BonsaiTree` with the new coordinates.

## Examples

``` r
# the same tree as a dendrogram
sc <- demo_single_cells(
  prepped = FALSE,
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 50L)
)
tree <- bonsai_sc(sc, .verbose = FALSE)
tree <- relayout_bonsai(tree, layout = "dendrogram")
plot(tree)


unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
