# Calculate recall at k and distance ratio

Helper function to compare the results of two
`SingleCellNearestNeighbour` against each other. The first one can serve
as a reference (ground truth) and you can compare against the second
one.

## Usage

``` r
calc_knn_metrics(ref_knn, query_knn)
```

## Arguments

- ref_knn:

  The reference `SingleCellNearestNeighbour`.

- query_knn:

  The query `SingleCellNearestNeighbour`.

## Value

A list with:

- matches - The intersecting indices between the reference and query kNN
  for each sample. In an ideal match up should be equal to k.

- distance_ratio - The distance ratio. Calculates
  `sum(dist_query) / sum(dist_ref)` per sample. Indicates how much worse
  the reference is.

- final_recall - The final recall across all samples.

- final_ratio - The final distance ratio across all samples.
