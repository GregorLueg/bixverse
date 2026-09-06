# Detect likely isotype-control features by name pattern

Detect likely isotype-control features by name pattern

## Usage

``` r
detect_adt_isotypes(feature_names, pattern = "isotype")
```

## Arguments

- feature_names:

  Character. ADT feature names (matrix colnames).

- pattern:

  String. Case-insensitive regex. Defaults to `"isotype"`.

## Value

Character vector of matching names, for inspection before passing to
[`add_adt_counts_sc()`](https://gregorlueg.github.io/bixverse/reference/add_adt_counts_sc.md)
as `isotype_names`.

## Examples

``` r
# find the isotype controls among the ADT features
features <- c("CD3", "CD19", "IgG1_isotype", "IgG2a_isotype")
detect_adt_isotypes(features)
#> [1] "IgG1_isotype"  "IgG2a_isotype"
```
