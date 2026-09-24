# Calculate LISI scores (iLISI or cLISI)

Computes the Local Inverse Simpson's Index (LISI) on the kNN graph: the
effective number of labels in each cell's neighbourhood. On batch labels
this is iLISI, where higher means better mixing. On cell type labels it
is cLISI, where lower means cell types stay apart. Unlike kBET, LISI
does not compare against global proportions, so it also works on
graph-based corrections like BBKNN.

The normalised score follows scIB and lands in `[0, 1]`, higher is
better for both: iLISI as `(median - 1) / (n - 1)`, cLISI as
`(n - median) / (n - 1)`.

## Usage

``` r
calculate_lisi_sc(
  object,
  label_column,
  type = c("batch", "cell_type"),
  weighted = FALSE,
  perplexity = 30,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- label_column:

  String. The column with the batch or cell type labels in the obs data
  of the class.

- type:

  String. One of `c("batch", "cell_type")`. Decides which normalised
  score is reported. Defaults to `"batch"`.

- weighted:

  Boolean. Weight the neighbours with a perplexity-calibrated Gaussian
  kernel on the kNN distances, as in Korsunsky et al. If `FALSE`, all
  neighbours count equally. Defaults to `FALSE`.

- perplexity:

  Numeric. Perplexity for the weighted version. Defaults to `30`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `LisiScores` object with the following elements

- per_cell - Per-cell LISI scores in `[1, n_labels]`.

- mean_lisi - Mean LISI across all cells.

- median_lisi - Median LISI across all cells.

- lisi_norm - The normalised score in `[0, 1]`, higher is better.

- n_labels - Number of distinct labels.

- type - `"batch"` (iLISI) or `"cell_type"` (cLISI).

## References

Korsunsky, et al., Nat. Methods, 2019; Luecken, et al., Nat. Methods,
2022

## Examples

``` r
# iLISI over the kNN graph
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_lisi_sc(
  sc,
  label_column = "batch_index",
  .verbose = FALSE
)
#> iLISI (batch)
#>   Cells: 600 | Labels: 3
#>   Mean LISI:    1.8192
#>   Median LISI:  1.8000
#>   Normalised:   0.4000 (0 = worst, 1 = best)

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
