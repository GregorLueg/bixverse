# Get factor matrices from a BulkModuleResult

Returns one factor matrix by key, or the whole named list of factor
matrices if `which` is `NULL`. Keys are method-specific; see the
`factors` argument of
[`new_bulk_module_result()`](https://gregorlueg.github.io/bixverse/reference/new_bulk_module_result.md)
for the shared vocabulary.

## Usage

``` r
get_factors(object, which = NULL)
```

## Arguments

- object:

  A `BulkModuleResult`.

- which:

  String or `NULL`. Name of the factor matrix to return. If `NULL`,
  returns the whole list.

## Value

A matrix, the named list of matrices, or `NULL` (with warning) if
`which` is not among the stored factor keys.
