# Get the final results from the class

Get the final results from `BixverseBaseClass` class (or child classes).

## Usage

``` r
get_results(object)

get_results.DialogueResult(object, ...)
```

## Arguments

- object:

  The underlying
  [`BixverseBaseClass()`](https://gregorlueg.github.io/bixverse/reference/BixverseBaseClass.md)
  class. The class functionality is usually inherited by other S7
  classes in `bixverse`.

- ...:

  Unused, present so the S3 methods sharing this page match the generic.

## Value

Returns the final results if any have been stored in the class.
