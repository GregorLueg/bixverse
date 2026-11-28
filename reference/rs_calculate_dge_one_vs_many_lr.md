# Calculate one-vs-many specific markers with logistic regression p-values

**\[experimental\]** The logistic regression twin of
[`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md):
AUROC summaries as there, p-values from the likelihood ratio test.

## Usage

``` r
rs_calculate_dge_one_vs_many_lr(
  f_path_cell,
  f_path_gene,
  cell_groups,
  references,
  min_prop,
  alternative,
  verbose
)
```

## Arguments

- f_path_cell:

  String. Path to the `counts_cells.bin` file, for the library sizes.

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- cell_groups:

  List. List of integer vectors, each containing the index positions
  (0-indexed) of the cells of one group.

- references:

  Integer. Index positions (0-indexed) into `cell_groups` of the
  reference groups to report.

- min_prop:

  Numeric. Minimum proportion of expression in at least one of the
  groups to be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

The list described in
[`rs_calculate_dge_one_vs_many()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_many.md).
