# Resolve a modality against the object it is used on

On a `SingleCellsAtac` `"atac"` names the object's own cache, i.e.
`"rna"`. Any modality other than `"rna"` on an object without other
modalities is an error.

## Usage

``` r
.resolve_modality(x, modality)
```

## Arguments

- x:

  The single cell object.

- modality:

  String. The requested modality.

## Value

String. The modality to use on `x`.
