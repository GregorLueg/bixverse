# Assert per-cell ScType parameters

Checkmate extension for asserting the per-cell ScType parameters as
returned by
[`params_sctype_cells()`](https://gregorlueg.github.io/bixverse/reference/params_sctype_cells.md).

## Usage

``` r
assertSctypeCellParams(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The list to assert.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns `x` if the assertion is successful.
