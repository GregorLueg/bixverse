# Calculate batch LISI scores

Computes the Local Inverse Simpson's Index (LISI) on batch labels using
the kNN graph. LISI measures the effective number of batches represented
in each cell's neighbourhood. Under perfect mixing, LISI equals the
number of batches. Under no mixing, LISI equals 1. Unlike kBET, LISI
does not compare against global batch proportions, making it suitable
for graph-based correction methods like BBKNN.

## Usage

``` r
calculate_batch_lisi_sc(object, batch_column, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `BatchLisiScores` object with the following elements

- per_cell - Per-cell LISI scores in `[1, n_batches]`.

- mean_lisi - Mean LISI across all cells.

- median_lisi - Median LISI across all cells.

- n_batches - Number of batches in the data.

## References

Korsunsky, et al., Nat. Methods, 2019

## Examples

``` r
# batch LISI over the kNN graph
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_batch_lisi_sc(
  sc,
  batch_column = "batch_index",
  .verbose = FALSE
)
#> Batch LISI Scores
#>   Cells: 600 | Batches: 3
#>   Mean LISI:    1.8192 (1 = no mixing, 3 = perfect mixing)
#>   Median LISI:  1.8000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
