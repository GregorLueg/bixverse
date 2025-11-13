# Create sparse dgRMatrix matrices for raw and norm counts

Create sparse dgRMatrix matrices for raw and norm counts

## Usage

``` r
get_meta_cell_matrices(meta_cell_data)
```

## Arguments

- meta_cell_data:

  Named list. This contains the indptr, indices and data for both raw
  counts and norm counts.

## Value

A list of two items

- raw - Sparse matrix representing the raw counts.

- norm - Sparse matrix representing the norm counts.
