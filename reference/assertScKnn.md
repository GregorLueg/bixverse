# Assert single cell kNN parameters

Checkmate extension for asserting the single cell kNN parameters.

## Usage

``` r
assertScKnn(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The list to check/assert.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertionFunction()`](https://mllg.github.io/checkmate/reference/makeAssertion.html).

## Value

Invisibly returns the checked object if the assertion is successful.
