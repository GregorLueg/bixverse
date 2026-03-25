# Calculate kBET scores

This function calculates the k-nearest neighbour batch-effect test
(kBET). Briefly, the function leverages a Chi-Square statistic to
calculate the differences in batch proportions observed in the
neighbourhood of a given cell with the overall batch proportions. If the
test is significant for that cell it indicates poor mixing for that cell
specifically. Large number of positive tests indicate bad mixing
overall. For more details, please see Büttner et al.

## Usage

``` r
calculate_kbet_sc(object, batch_column, threshold = 0.05, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- threshold:

  Numeric. Number between 0 and 1. Below this threshold, the test is
  considered significant. Defaults to `0.05`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A `KbetScores` object with the following elements

- kbet_score - Proportion of significant tests over all cells. 0
  indicates perfect mixing, 1 indicates no mixing between batches.

- significant_tests - Logical vector indicating for which cells the test
  was below the threshold.

- p_values - The p-values from the Chi-Square test.

- chi_square_stats - Per-cell Chi-Square statistics.

- mean_chi_square - Mean Chi-Square statistic across all cells.

- median_chi_square - Median Chi-Square statistic across all cells.

- threshold - The significance threshold used.

- n_batches - Number of batches in the data.

## References

Büttner, et al., Nat. Methods, 2019
