# Helper to compare kNN graphs

**\[experimental\]** Compare two kNN graphs and return the distance
ratios and overlaps of k-nearest neighbours between them.

## Usage

``` r
rs_compare_knn(knn_data_a, knn_data_b)
```

## Arguments

- knn_data_a:

  Named list. This contains the kNN data (including distances) of the
  first kNN graph. This one will be treated as the ground truth

- knn_data_b:

  Named list. This contains the kNN data (including distances) of the
  second kNN graph.

## Value

A list with the following elements:

- all_matches - Integer vector. Number of shared neighbours per sample.

- all_ratios - Numerical vector. Ratio of the summed distances (b / a)
  per sample. Samples whose summed distance in a is ~0 are skipped, so
  this can be shorter than `all_matches`.

- final_recall - The mean recall across all samples, with a as the
  ground truth.

- final_ratio - The mean distance ratio across the retained samples.
