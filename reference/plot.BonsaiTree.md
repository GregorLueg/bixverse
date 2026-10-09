# Plot a Bonsai tree

Draws every branch of the tree and a point per leaf, a cell or a
metacell. Radial layouts are drawn with straight branches at a fixed
aspect ratio, the dendrogram with right angles.

## Usage

``` r
# S3 method for class 'BonsaiTree'
plot(
  x,
  colour_by = NULL,
  layout = NULL,
  hyperbolic = NULL,
  point_size = 0.5,
  size_by_cells = FALSE,
  edge_colour = "grey60",
  ...
)
```

## Arguments

- x:

  A `BonsaiTree` object.

- colour_by:

  Optional vector with one value per cell, in the cell order of the tree
  (`x$nodes$cell_id[x$nodes$is_leaf]`). Numeric values get a continuous
  scale, anything else a discrete one.

- layout:

  Optional string. One of
  `c("equal_angle", "equal_daylight", "dendrogram")`. If it differs from
  the stored layout, the tree is laid out again first, see
  [`relayout_bonsai()`](https://gregorlueg.github.io/bixverse/reference/relayout_bonsai.md).

- hyperbolic:

  Optional boolean. As `layout`, for the hyperbolic projection.

- point_size:

  Numeric. Size of the leaf points; the smallest size if
  `size_by_cells = TRUE`.

- size_by_cells:

  Boolean. Scale each leaf point by the number of cells behind it, which
  only differs between leaves for metacells.

- edge_colour:

  String. Colour of the branches.

- ...:

  Additional arguments (unused; required by the S3 generic).

## Value

A `ggplot2` object.

## Examples

``` r
# the tree coloured by cell type
sc <- demo_single_cells(
  prepped = FALSE,
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 50L)
)
tree <- bonsai_sc(sc, .verbose = FALSE)
plot(tree, colour_by = get_sc_obs(sc, filtered = TRUE)$cell_grp)


unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
