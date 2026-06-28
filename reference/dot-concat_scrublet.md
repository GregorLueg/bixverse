# Concatenate per-group Scrublet results into a single `ScrubletRes` object

Cell-level vectors are re-ordered by the original cell index. Per-group
scalar summaries (threshold, rates) are retained as named vectors.

## Usage

``` r
.concat_scrublet(group_results, group_by, return_combined_pca, return_pairs)
```

## Arguments

- group_results:

  Named list of `ScrubletRes` objects.

- group_by:

  Character. Name of the grouping column; stored as an attribute on the
  returned object.

- return_combined_pca:

  Logical. If `TRUE`, PCA results are included per group.

- return_pairs:

  Logical. If `TRUE`, simulated doublet pair indices are included per
  group.

## Value

A `ScrubletRes` object with attributes `cell_indices`, `grouped`, and
`group_by_col`.
