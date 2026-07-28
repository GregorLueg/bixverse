# Apply a pipeline independently to each group of a `SingleCells` object

Splits `object` by `group_col`, applies `pipeline` to each subset, and
returns a named list of processed `SingleCellsSubset`s. Useful for
per-sample / per-cell-type re-analysis where the same chain (HVG, PCA,
neighbours, clusters, ...) is run on each group, e.g. sample-pure
metacell generation followed by an external merge.

## Usage

``` r
apply_pipeline_per_group(pipeline, object, group_col, groups = NULL)
```

## Arguments

- pipeline:

  `ScPipeline`.

- object:

  `SingleCells`.

- group_col:

  String. Column in obs used to split.

- groups:

  Optional character vector. Restrict to these group values; if `NULL`,
  all unique values of `group_col` are used.

## Value

Named list of `SingleCellsSubset` objects, names being the group values.
