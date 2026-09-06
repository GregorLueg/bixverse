# Helper function to write data to h5ad format

This is a helper to write synthetic data to h5ad file. This version will
write the data into the common compressed sparse data format.

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

## Examples

``` r
# round trip synthetic counts through a sparse h5ad
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
get_h5ad_dimensions(f_path)$dims
#> obs var 
#> 200  40 

unlink(f_path)
```
