# Get the parameters that were used.

Extracts parameters from the `bixverse_base_class` class (or child
classes) and has options to return (pretty) JSONs.

## Usage

``` r
get_params(object, to_json = FALSE, pretty_json = FALSE)
```

## Arguments

- object:

  A class within bixverse that inherits from
  [`bixverse_base_class()`](bixverse_base_class.md).

- to_json:

  Shall the params be returned as a JSON string.

- pretty_json:

  Shall the params be returned as a pretty JSON string.

## Value

Depending on parameters either the R list or a (pretty) JSON string.
