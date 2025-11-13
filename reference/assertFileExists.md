# Assert that files exist

Checkmate extension for asserting if files exist in the directory.

## Usage

``` r
assertFileExists(x, file_names, .var.name = checkmate::vname(x), add = NULL)
```

## Arguments

- x:

  String. The directory to check the files for.

- file_names:

  String. Vector of names of the expected files in this directory.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly returns the checked object if the assertion is successful.
