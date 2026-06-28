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

- all_matches - Matching neighbours for this sample.

- all_ratios - Distance ratio for this sample (with b / a).

- final_recall - The final recall of assuming a being the ground truth
  across all samples

- final_ratio - The final distance ratio across all samples
