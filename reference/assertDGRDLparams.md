# Assert DGRDL params

Assert DGRDL params

## Usage

``` r
assertDGRDLParams(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The object to check.

- .var.name:

  Name of the checked object to print in assertions.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns the checked object if the assertion is successful.
