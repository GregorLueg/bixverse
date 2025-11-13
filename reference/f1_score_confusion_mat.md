# F1 scores on top of a confusion matrix

Helper function to check for expected clustering vs actual clustering.

## Usage

``` r
f1_score_confusion_mat(clusters_a, clusters_b)
```

## Arguments

- clusters_a:

  String or factor. The clustering of algorithm 1.

- clusters_b:

  String or factor. The clustering of algorithm 2.

## Value

Named vector with the F1 scores between the two clustering algorithms.
