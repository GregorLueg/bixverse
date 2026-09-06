# Remove the KNN data

Remove the KNN data

## Usage

``` r
remove_knn(x, ...)

## S7 method for class <bixverse::MetaCells>
remove_knn(x, ...)

# S3 method for class 'ScCache'
remove_knn(x, ...)

## S7 method for class <bixverse::SingleCells>
remove_knn(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
remove_knn(x, ...)
```

## Arguments

- x:

  An object from which to remove the kNN data.

- ...:

  Other parameters.

## Examples

``` r
# drop the cached kNN, for example before rebuilding it
sc <- demo_single_cells()
sc <- remove_knn(sc)
is.null(get_knn_obj(sc))
#> [1] TRUE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
