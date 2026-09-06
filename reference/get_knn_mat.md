# Get the KNN matrix

Getter for an integer matrix of samples x neighbours.

## Usage

``` r
get_knn_mat(x, ...)

## S7 method for class <bixverse::MetaCells>
get_knn_mat(x, ...)

# S3 method for class 'ScCache'
get_knn_mat(x, ...)

## S7 method for class <bixverse::SingleCells>
get_knn_mat(x, ...)

# S3 method for class 'SingleCellNearestNeighbour'
get_knn_mat(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_knn_mat(x, ...)
```

## Arguments

- x:

  An object to get the kNN matrix from.

- ...:

  Other parameters.

## Examples

``` r
# cells x neighbours
sc <- demo_single_cells()
dim(get_knn_mat(sc))
#> [1] 500  15

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
