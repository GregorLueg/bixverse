# Helper function to generate manifoldsR nearest neighbours

Helper function to generate manifoldsR nearest neighbours

## Usage

``` r
.get_manifoldsr_knn(x, modality = c("rna", "adt", "atac"))
```

## Arguments

- x:

  `SingleCells`, `MetaCells` object from which to extract the kNN data.

- modality:

  String. One of `c("rna", "adt", "atac")`.

## Value

The manifoldsR `NearestNeigbours` if possible
