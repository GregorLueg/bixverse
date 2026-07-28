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

  `MetaCells` class.

- original_cell_type:

  Character vector. The original cell type annotations, in the row order
  of the full (unfiltered) obs table of the object the meta cells came
  from, i.e. `get_sc_obs(x)$<column>`. Meta cell memberships are stored
  as 1-indexed positions in that space, so a vector restricted to the
  QC-passing cells will silently give wrong purities.

## Value

The `MetaCells` with an added columns to the observation table with the
purity measures
