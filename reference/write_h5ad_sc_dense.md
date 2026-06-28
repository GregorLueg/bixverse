# Helper function to write data to a dense h5ad file

Helper to write synthetic data to h5ad with `/X` stored as a dense 2D
dataset (cells x genes). Useful for exercising the dense ingestion path.

## Usage

``` r
write_h5ad_sc_dense(
  f_path,
  counts,
  obs,
  var,
  overwrite = TRUE,
  .verbose = TRUE
)
```

## Arguments

- f_path:

  String. The filepath to which to save the data

- counts:

  Matrix or sparse matrix; sparse input is densified.

- obs:

  data.table. The observations. Needs `nrow(obs) == nrow(counts)`.

- var:

  data.table. The variable data. Needs `nrow(var) == ncol(counts)`.

- overwrite:

  Boolean. Shall any found h5ad file be overwritten.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Returns invisible
