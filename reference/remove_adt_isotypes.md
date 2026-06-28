# Return the ADT feature names removing the isotypes

Return the ADT feature names removing the isotypes

## Usage

``` r
remove_adt_isotypes(feature_names, pattern = "isotype")
```

## Arguments

- feature_names:

  Character. ADT feature names (matrix colnames or from
  [`get_adt_names()`](https://gregorlueg.github.io/bixverse/reference/get_adt_names.md)).

- pattern:

  String. Case-insensitive regex. Defaults to `"isotype"`.

## Value

Character vector of ADT features, but anything with `"isotype"`.
