# Transform Seurat raw counts into a List

Transform Seurat raw counts into a List

## Usage

``` r
get_seurat_counts_to_list(seurat_obj)
```

## Arguments

- seurat_obj:

  `Seurat` class. The class from which to extract the counts from and
  transform to a list.

## Value

A list with the following elements

- indptr - Index pointers of the sparse data.

- indices - Indices of the data.

- data - The underlying data.

- format - String that defines if the data is CSR or CSC.

- nrow - The number of rows.

- ncol - The number of columns.
