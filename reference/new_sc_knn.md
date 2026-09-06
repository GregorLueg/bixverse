# Helper function to generate kNN data with distances

Wrapper class that stores kNN data for subsequent usage in various other
functions, methods and helpers.

## Usage

``` r
new_sc_knn(knn_data, used_cells)
```

## Arguments

- knn_data:

  Named list with the following items:

  - indices - Integer matrix containing the indices of the nearest
    neighbours (0-indexed).

  - dist - Numerical matrix containing the distances to the nearest
    neighbours.

  - dist_metric - String. Distance metric used.

- used_cells:

  Character vector. The cells used to generate the kNN graph with the
  distances.

## Value

Generates the `SingleCellNearestNeighbour` class.

## Examples

``` r
# rebuild the kNN wrapper from an existing graph
sc <- demo_single_cells()
knn <- generate_knn_sc(sc, .validate_index = FALSE, .verbose = FALSE)
new_sc_knn(
  knn_data = list(
    indices = get_knn_mat(knn),
    dist = get_knn_dist(knn),
    dist_metric = "euclidean"
  ),
  used_cells = get_cell_names(sc)
)
#> SingleCellNearestNeighbour: 500 cells, k = 15
#>   Distance metric: euclidean
#>   Index range: [0, 499]
#>   Distance range: [1.1849, 4.2876]

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
