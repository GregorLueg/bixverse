# Get the PCA singular values

Returns the PCA singular values (can be useful to assess cumulative
variance explained). This function is used for the single cell-related
classes and methods.

## Usage

``` r
get_pca_singular_val(x, ...)

## S7 method for class <bixverse::MetaCells>
get_pca_singular_val(x, ...)

# S3 method for class 'ScCache'
get_pca_singular_val(x, ...)

## S7 method for class <bixverse::SingleCells>
get_pca_singular_val(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_pca_singular_val(x, ...)
```

## Arguments

- x:

  An object to get PCA singular values from.

- ...:

  Other parameters.

## Value

The PCA singular values from the object (if found).

## Examples

``` r
# singular values, largest first
sc <- demo_single_cells()
head(get_pca_singular_val(sc), 3)
#> [1] 47.75282 44.16640 25.33718

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
