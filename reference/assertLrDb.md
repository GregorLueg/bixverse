# Assert a ligand-receptor database

Checkmate extension for asserting a ligand-receptor database as returned
by
[`get_cellphonedb_db()`](https://gregorlueg.github.io/bixverse/reference/get_cellphonedb_db.md).

## Usage

``` r
assertLrDb(x, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  The data.table to assert.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns `x` if the assertion is successful.
