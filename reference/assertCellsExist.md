# Assert neighbour generation parameters

Checkmate extension for asserting if the prodivided cell names

## Usage

``` r
assertCellsExist(x, cell_names, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The `single_cell_exp` object to check/assert.

- cell_names:

  String. The provided cell names.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns the checked object if the assertion is successful.
