# Assert cell marker list

Checkmate extension for asserting a cell marker list as returned by
[`prepare_cell_markers()`](https://gregorlueg.github.io/bixverse/reference/prepare_cell_markers.md).

## Usage

``` r
assertCellMarkerList(x, .var.name = checkmate::vname(x), add = NULL)
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
