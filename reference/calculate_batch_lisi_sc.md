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

  `SingleCells` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A `BatchLisiScores` object with the following elements

- per_cell - Per-cell LISI scores in `[1, n_batches]`.

- mean_lisi - Mean LISI across all cells.

- median_lisi - Median LISI across all cells.

- n_batches - Number of batches in the data.

## References

Korsunsky, et al., Nat. Methods, 2019
