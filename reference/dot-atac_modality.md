# The modality name of the ATAC cache of an object

The modality name of the ATAC cache of an object

## Usage

``` r
.atac_modality(x)
```

## Arguments

- x:

  `SingleCellsAtac` or `SingleCellsMultiModal`.

## Value

String. `"rna"` for `SingleCellsAtac`, whose own cache is the ATAC one,
`"atac"` otherwise.
