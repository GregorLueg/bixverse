# Helper function to write data to h5ad format

This is a helper to write synthetic data to h5ad file.

## Usage

``` r
write_h5ad_sc(f_path, counts, obs, var, overwrite = TRUE, .verbose = TRUE)
```

## Arguments

- f_path:

  String. The filepath to which to save the data

- counts:

  Sparse matrix. Needs to be of class `dgRMatrix` or `dgCMatrix`.

- obs:

  data.table. The observations. Needs to have
  `nrow(obs) == nrow(counts)`.

- var:

  data.table. The variable data. Needs to have
  `ncol(var) == ncol(counts)`.

- overwrite:

  Boolean. Shall any found h5ad file be overwritten.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Returns invisible
