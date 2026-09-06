# Add additional embeddings to the class

Add additional embeddings to the class

## Usage

``` r
set_embedding(x, embd, name, ...)

## S7 method for class <bixverse::MetaCells>
set_embedding(x, embd, name, ...)

# S3 method for class 'ScCache'
set_embedding(x, embd, name, ...)

## S7 method for class <bixverse::SingleCells>
set_embedding(x, embd, name, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_embedding(x, embd, name, ...)
```

## Arguments

- x:

  An object to add the singular values for.

- embd:

  Numerical matrix representing the additional embedding.

- name:

  String. Name of the embedding.

- ...:

  Other parameters.

## Examples

``` r
# the first two PCs stored as an embedding in their own right
sc <- demo_single_cells()
sc <- set_embedding(sc, get_pca_factors(sc)[, 1:2], name = "pca_2d")
get_available_embeddings(sc)
#> [1] "pca"    "pca_2d"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
