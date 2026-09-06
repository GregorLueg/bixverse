# Set the HVG genes

Stores within the class the index positions of the HVG. This is used for
the single cell-related classes and methods.

## Usage

``` r
set_hvg(x, hvg)

## S7 method for class <bixverse::MetaCells>
set_hvg(x, hvg)

# S3 method for class 'ScMap'
set_hvg(x, hvg)

## S7 method for class <bixverse::SingleCells>
set_hvg(x, hvg)

## S7 method for class <bixverse::SingleCellsSubset>
set_hvg(x, hvg)
```

## Arguments

- x:

  An object to set the HVGs for

- hvg:

  String or integer. The names or indices of the highly variable genes.

## Examples

``` r
# HVGs picked by hand instead of by find_hvg_sc()
sc <- demo_single_cells(prepped = FALSE)
sc <- set_hvg(sc, get_gene_names(sc)[1:20])
head(get_hvg(sc))
#> [1] 0 1 2 3 4 5

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
