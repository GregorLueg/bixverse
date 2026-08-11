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

  Character vector. The original cell type annotations of the object the
  meta cells came from. Either in the row order of its full (unfiltered)
  obs table, i.e. `get_sc_obs(x)$<column>`, or of the QC-passing cells
  only, i.e. `get_sc_obs(x, filtered = TRUE)$<column>`. Which one you
  passed is inferred from the length, so a vector matching neither is an
  error rather than a silently wrong purity.

## Value

The `MetaCells` with an added columns to the observation table with the
purity measures
