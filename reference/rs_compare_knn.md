# Helper to compare kNN graphs

Helper to compare kNN graphs

## Usage

``` r
rs_compare_knn(knn_data_a, knn_data_b)
```

## Arguments

- knn_mat_a:

  Integer matrix. The indices of the first kNN graph to compare. Should
  be samples x neighbours. This will be treated as ground truth.

- knn_mat_b:

  Integer matrix. The indices of the second kNN graph to compare. Should
  be samples x neighbours.

- knn_dist_a:

  Numeric matrix.

## Value

A list with the following elements:

- all_matches - Matching neighbours for this sample.

- all_ratios - Distance ratio for this sample (with b / a).

- final_recall - The final recall of assuming a being the ground truth
  across all samples

- final_ratio - The final distance ratio across all samples
