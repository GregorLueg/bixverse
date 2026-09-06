# Helper function to get the dimensions and storage format

Distinguishes sparse (`CSR`/`CSC`) from dense storage. For dense `X` the
orientation is inferred by comparing the matrix dims against the obs and
var lengths (`DENSE_ROW` = cells x genes, `DENSE_COL` = genes x cells).
Ties (`no_obs == no_var`) fall back to the AnnData convention
(`DENSE_ROW`).

## Usage

``` r
get_h5ad_dimensions(f_path)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

## Value

A list with `dims` (named integer `c(obs, var)`) and `type` (one of
`"CSR"`, `"CSC"`, `"DENSE_ROW"`, `"DENSE_COL"`).

## Examples

``` r
# dimensions and storage layout without reading the counts
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
get_h5ad_dimensions(f_path)
#> $dims
#> obs var 
#> 200  40 
#> 
#> $type
#> [1] "CSR"
#> 

unlink(f_path)
```
