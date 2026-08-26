# Assert LDA parameters

Checkmate extension for asserting the LDA solver parameters, see
[`params_lda()`](https://gregorlueg.github.io/bixverse/reference/params_lda.md).

## Usage

``` r
assertLdaParams(x, .var.name = checkmate::vname(x), add = NULL)
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
