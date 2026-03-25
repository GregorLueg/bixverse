# Return the resolution results

Getter function to get the resolution results (if available).

## Usage

``` r
get_resolution_res(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

If resolution results were found, returns the data.table. Otherwise,
throws a warning and returns NULL.
