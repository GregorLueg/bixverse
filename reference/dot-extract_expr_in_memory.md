# Extract dense expression for in-memory count matrices

Extract dense expression for in-memory count matrices

## Usage

``` r
.extract_expr_in_memory(mat, features, scale = FALSE, clip = NULL)
```

## Arguments

- mat:

  Numeric matrix. Cells x features, with cell ids as row names and
  feature ids as column names.

- features:

  Character vector. Feature ids to extract (pre-matched).

- scale:

  Boolean. Whether to z-score each feature across cells.

- clip:

  Optional numeric. Clip the z-scores if `scale = TRUE`.

## Value

A data.table with a `cell_id` column and one column per feature.
