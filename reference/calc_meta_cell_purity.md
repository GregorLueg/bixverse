# Calculate meta cell purity

A potential metric to see how well the meta cells are aggregated is
their cell type purity. This helper function helps to plot the meta-cell
purity based on annotated cell types. These can be also just
unsupervised memberships to graph-based clustering, etc.

## Usage

``` r
calc_meta_cell_purity(object, original_cell_type)
```

## Arguments

- object:

  `meta_cells` class.

- original_cell_type:

  Character vector. The original cell type annotations. The indices need
  to match with the original cell indices used to generate the
  meta-cells! (1-indexed)

## Value

The `meta_cells` with an added columns to the observation table with the
purity measures
