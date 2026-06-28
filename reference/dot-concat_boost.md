# Concatenate per-group boost results into a single `BoostRes` object

Cell-level vectors are re-ordered by the original cell index.

## Usage

``` r
.concat_boost(group_results, group_by)
```

## Arguments

- group_results:

  Named list of `BoostRes` objects.

- group_by:

  Character. Name of the grouping column; stored as an attribute on the
  returned object.

## Value

A `BoostRes` object with attributes `cell_indices`, `grouped`, and
`group_by_col`.
