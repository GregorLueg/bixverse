# Get k-means cluster assignments from a fast cluster result

Get k-means cluster assignments from a fast cluster result

## Usage

``` r
get_kmeans_clusters(x)

# S3 method for class 'SingleCellFastClusters'
get_kmeans_clusters(x)
```

## Arguments

- x:

  `SingleCellFastClusters` object.

## Examples

``` r
# the centroid each cell was assigned to
sc <- demo_single_cells()
res <- fast_cluster_sc(
  sc,
  resolutions = 1.0,
  n_centroids = 30L,
  return_kmeans = TRUE,
  .verbose = FALSE
)
head(get_kmeans_clusters(res))
#> [1]  5  9 20 19  2 20

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
