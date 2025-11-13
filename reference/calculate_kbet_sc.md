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

  `single_cell_exp` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- threshold:

  Numeric. Number between 0 and 1. Below this threshold, the test is
  considered significant. Defaults to `0.05`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A list with the following elements

- kbet_score - Number of significant tests over all cells. 0 indicates
  perfect mixing, 1 indicates basically zero mixing between batches.

- significant_tests - Boolean indicating for which cells the statistic
  was below the threshold

- chisquare_pvals - The p-values of the ChiSquare test.

## References

Büttner, et al., Nat. Methods, 2019
