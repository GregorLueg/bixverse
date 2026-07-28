# bixverse single cell subset class

Subset view onto a
[`SingleCells()`](https://gregorlueg.github.io/bixverse/reference/SingleCells.md)
object, restricted to cells belonging to a single level of a grouping
variable. The Rust count connection is shared with the parent (no data
copy). `obs_table` and `var_table` are held in memory; `sc_map` is
rebuilt to point only at the subset cells but stays in the original
index space so Rust calls remain valid without further translation.

## Usage

``` r
SingleCellsSubset(sc_object, grouping_column, group)
```

## Value

A `SingleCellsSubset` object.

## Properties

- count_connection:

  Shared Rust pointer to the on-disk counts.

- dir_data:

  Directory holding the binary count files.

- obs_table:

  Subset obs (rows for the chosen group only). `cell_idx` keeps the
  original 1-indexed position in the parent.

- var_table:

  Variable/feature table (unchanged from parent).

- grouping_column:

  Column in obs used to define the subset.

- group:

  Value of `grouping_column` represented by this subset.

- sc_cache:

  Fresh `ScCache` for subset-specific PCA, kNN, sNN, embeddings.

- sc_map:

  `ScMap` restricted to the subset cells. `cell_mapping` stays 1-indexed
  and `cells_to_keep_idx` stays 0-indexed, both in the original parent
  index space.

- subset_to_original:

  Integer vector. 1-indexed original cell positions, in subset row
  order. `subset_to_original[i]` is the parent position of subset row
  `i`.

- dims:

  `c(n_cells_subset, n_genes)`.
