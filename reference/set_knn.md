# Set/add KNN

Set/add KNN

## Usage

``` r
set_knn(x, knn, ...)

## S7 method for class <bixverse::MetaCells>
set_knn(x, knn, ...)

# S3 method for class 'ScCache'
set_knn(x, knn, ...)

## S7 method for class <bixverse::SingleCells>
set_knn(x, knn, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_knn(x, knn, ...)
```

## Arguments

- x:

  An object to add the KNN data to

- knn:

  `SingleCellNearestNeighbour` class to add to the classes.

- ...:

  Other parameters.

## Examples

``` r
# a kNN built outside the object put into the cache
sc <- demo_single_cells()
knn <- generate_knn_sc(sc, .validate_index = FALSE, .verbose = FALSE)
sc <- set_knn(sc, knn)
dim(get_knn_mat(sc))
#> [1] 500  15

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
