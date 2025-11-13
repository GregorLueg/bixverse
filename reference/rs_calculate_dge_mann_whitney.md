# Calculate DGEs between cells based on Mann Whitney stats

The function will take two sets of cell indices and calculate the
differential gene expression based on the Mann Whitney test between the
two groups.

## Usage

``` r
rs_calculate_dge_mann_whitney(
  f_path,
  cell_indices_1,
  cell_indices_2,
  min_prop,
  alternative,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- cell_indices_1:

  Integer. Index positions (0-indexed) of the cells of group 1.

- cell_indices_2:

  Integer. Index positions (0-indexed) of the cells of group 2.

- min_prop:

  Minimum proportion of expression in at least one of the two groups to
  be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements

- lfc - Log fold changes between the two groups.

- prop1 - Proportion of cells expressing the gene in group 1.

- prop2 - Proportion of cells expressing the gene in group 2.

- z_scores - Z-scores based on the Mann Whitney statistic.

- p_values - P-values of the Mann Whitney statistic.

- fdr - False discovery rate after BH adjustment

- genes_to_keep - Boolean indicating which genes were tested.
