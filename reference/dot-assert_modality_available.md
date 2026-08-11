# Guard against non-RNA modalities on single-modality classes

Guard against non-RNA modalities on single-modality classes

## Usage

``` r
.assert_modality_available(object, modality)
```

## Arguments

- object:

  One of the single cell classes.

- modality:

  String. The requested modality.

## Value

`NULL`, invisibly. Errors if the modality is unavailable.
