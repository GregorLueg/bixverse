# Merge obs columns from subsets back into the parent object

Takes columns computed on one or more
[`SingleCellsSubset()`](https://gregorlueg.github.io/bixverse/reference/SingleCellsSubset.md)
objects and writes them into the parent's obs table in the DuckDB,
joining on `cell_idx`. The join is a left join, so every parent cell
that is not part of any of the provided subsets ends up as `NA`.

Pass a single subset or a list of them, e.g. the output of
[`apply_pipeline_per_group()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline_per_group.md).
All subsets are written in a single join, and the function refuses to
run if two subsets claim the same cell.

## Usage

``` r
merge_subset_obs(
  object,
  subsets,
  cols = NULL,
  new_names = NULL,
  prefix_values = FALSE,
  overwrite = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class. The parent the subsets were derived from.

- subsets:

  A `SingleCellsSubset` or a list of them. Every element must originate
  from `object`; this is verified against the cell names.

- cols:

  Optional character vector. Columns in the subset obs tables to merge.
  If `NULL`, every column that is present in all subsets but absent from
  the parent obs table is taken, which after a pipeline run is exactly
  the set of newly generated columns.

- new_names:

  Optional character vector. Names to give the merged columns in the
  parent obs table. Same length as `cols`, which must be given
  explicitly if this is used.

- prefix_values:

  Boolean. Prefix every merged value with the subset's group, i.e.
  `"<group>_<value>"`. Coerces the columns to character, so this only
  makes sense for discrete labels. Needed when merging a list of subsets
  into a shared column, since sub-cluster `1` of one group and
  sub-cluster `1` of another are otherwise indistinguishable.

- overwrite:

  Boolean. Allow merged columns to replace existing columns of the same
  name in the parent obs table. Defaults to `FALSE`, in which case a
  name clash is an error.

- .verbose:

  Boolean. Controls verbosity.

## Value

The parent `SingleCells` object with the merged columns added to the obs
table in the DuckDB.

## See also

[`add_sc_new_obs()`](https://gregorlueg.github.io/bixverse/reference/add_sc_new_obs.md)
for the equivalent on result objects that carry their own `cell_idx`,
such as
[`fast_cluster_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_cluster_sc.md).
