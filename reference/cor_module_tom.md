# Update the correlation matrix to a TOM

This function will update the correlation matrix to a topological
overlap matrix. It defaults to `"v2"` and the signed version, please see
[`calculate_tom()`](calculate_tom.md) for details.

## Usage

``` r
cor_module_tom(object, signed = TRUE, version = c("v2", "v1"), .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). You need to have
  applied [`cor_module_processing()`](cor_module_processing.md) before
  applying this function.

- signed:

  Boolean. Do you want to use the signed or unsigned version. Defaults
  to `TRUE`.

- version:

  String. One of `c("v2", "v1")`. Defaults to `"v2"`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.
