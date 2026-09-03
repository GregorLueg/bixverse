# Pull the raw counts out of a Seurat object

Pull the raw counts out of a Seurat object

## Usage

``` r
get_seurat_counts(seurat_obj)
```

## Arguments

- seurat_obj:

  `Seurat` class. The class to extract the counts from.

## Value

The raw counts as a `dgRMatrix` of cells x genes.
