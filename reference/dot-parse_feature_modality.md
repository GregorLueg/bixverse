# Parse an optional modality suffix from a feature id

Parse an optional modality suffix from a feature id

## Usage

``` r
.parse_feature_modality(feature, default_modality)
```

## Arguments

- feature:

  String. Feature id, optionally suffixed with `_rna` or `_adt`.

- default_modality:

  String. Modality used when no suffix is present.

## Value

A list with `id` and `modality`.
