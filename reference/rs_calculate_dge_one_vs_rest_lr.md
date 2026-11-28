# Calculate one-vs-rest logistic regression DGEs for every cell group

**\[experimental\]** The logistic regression twin of
[`rs_calculate_dge_one_vs_rest()`](https://gregorlueg.github.io/bixverse/reference/rs_calculate_dge_one_vs_rest.md).

## Usage

``` r
rs_calculate_dge_one_vs_rest_lr(
  f_path_cell,
  f_path_gene,
  cell_groups,
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

- min_prop:

  Numeric. Minimum proportion of expression in the group or in the rest
  to be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with `group`, `gene_idx` (both 0-indexed), `lfc`, `prop1`,
`prop2`, `z_scores`, `p_values` and `fdr`, one entry per tested gene and
group, flattened group-major.
