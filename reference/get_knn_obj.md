# Get the KNN object

Returns the `SingleCellNearestNeighbour` from the object. This function
is used for the single cell-related classes and methods.

## Usage

``` r
get_knn_obj(x, ...)

## S7 method for class <bixverse::MetaCells>
get_knn_obj(x, ...)

# S3 method for class 'ScCache'
get_knn_obj(x, ...)

## S7 method for class <bixverse::SingleCells>
get_knn_obj(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_knn_obj(x, ...)
```

## Arguments

- x:

  An object to get the KNN class from.

- ...:

  Other parameters.

## Value

The `SingleCellNearestNeighbour` object.

## Examples

``` r
# the cached kNN, cells x neighbours
sc <- demo_single_cells()
dim(get_knn_mat(get_knn_obj(sc)))
#> [1] 500  15

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
