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

## Examples

``` r
# agreement between two clusterings of the same six samples
clusters_a <- c("c1", "c1", "c2", "c2", "c3", "c3")
clusters_b <- c("x", "x", "y", "y", "z", "x")
f1_score_confusion_mat(clusters_a, clusters_b)
#>  c1  c2  c3 
#> 0.8 1.0 0.4 
```
