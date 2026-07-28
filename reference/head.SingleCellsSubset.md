# head Method for SingleCellsSubset object

head Method for SingleCellsSubset object

## Usage

``` r
## S7 method for class <bixverse::SingleCellsSubset>
head(x, n = 6L, ...)
```

## Arguments

- x:

  An object of class `SingleCellsSubset`.

- n:

  Integer. Number of rows to return. Defaults to `6L`.

- ...:

  Additional arguments (currently not used).

## Value

A data.table with the first `n` rows of the obs table.
