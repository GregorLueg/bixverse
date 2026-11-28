# The embedding Harmony corrects

The embedding Harmony corrects

## Usage

``` r
.harmony_input(object, embd_to_use, modality)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- embd_to_use:

  String. `"pca"` or another cached embedding.

- modality:

  String. The resolved modality.

## Value

The cells x dimensions matrix, or `NULL` with a warning when the
embedding is absent.
