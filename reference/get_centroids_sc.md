# Get k-means centroids from a fast cluster result

Get k-means centroids from a fast cluster result

## Usage

``` r
get_centroids_sc(x)

# S3 method for class 'SingleCellFastClusters'
get_centroids_sc(x)
```

## Arguments

- x:

  `SingleCellFastClusters` object.

## Examples

``` r
# the k-means centroids the graph clustering was built on
sc <- demo_single_cells()
res <- fast_cluster_sc(
  sc,
  resolutions = 1.0,
  n_centroids = 30L,
  return_kmeans = TRUE,
  .verbose = FALSE
)
dim(get_centroids_sc(res))
#> [1] 30 10

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
