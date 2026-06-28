# Concatenate per-group scDblFinder results into a single `ScDblFinderRes` object

Cell-level vectors are re-ordered by the original cell index. Cluster
labels are prefixed with the group name to avoid collisions across
groups. If features were requested, per-group feature matrices are
row-bound and similarly reordered.

## Usage

``` r
.concat_scdblfinder(group_results, group_by)
```

## Arguments

- group_results:

  Named list of `ScDblFinderRes` objects.

- group_by:

  Character. Name of the grouping column; stored as an attribute on the
  returned object.

## Value

A `ScDblFinderRes` object with attributes `cell_indices`, `grouped`, and
`group_by_col`.
