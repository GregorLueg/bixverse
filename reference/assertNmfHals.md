# Assert NMF HALS parameters

Checkmate extension for asserting the NMF HALS parameters for single
cell.

## Usage

``` r
assertNmfHals(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The list to check/assert

- .var.name:

  Name of the checked object to print in assertions.

- add:

  Collection to store assertion messages.

## Value

Invisibly returns the checked object if the assertion is successful.
