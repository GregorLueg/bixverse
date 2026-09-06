# Set/add the MAGIC imputed layer

Set/add the MAGIC imputed layer

## Usage

``` r
set_magic(x, magic, ...)

# S3 method for class 'ScCache'
set_magic(x, magic, ...)

## S7 method for class <bixverse::SingleCells>
set_magic(x, magic, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_magic(x, magic, ...)
```

## Arguments

- x:

  An object to add the imputed layer to.

- magic:

  `ScMagic` class with the imputed counts.

- ...:

  Other parameters.

## Value

The object with the imputed layer attached.

## Examples

``` r
# the imputed layer taken out and put back
sc <- demo_single_cells()
sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
magic <- get_magic(sc)
sc <- set_magic(remove_magic(sc), magic)
dim(get_magic(sc)$data)
#> [1] 500   5

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
