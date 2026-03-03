# Assert label propagation parameters

Checkmate extension for asserting label propagation parameters.

## Usage

``` r
assertLabelPropParams(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The list to check/assert.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns the checked object if the assertion is successful.
