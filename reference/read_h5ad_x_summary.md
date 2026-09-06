# Read summary statistics from the X slot of an h5ad file

Read summary statistics from the X slot of an h5ad file

## Usage

``` r
read_h5ad_x_summary(f_path, n_sample = 10000L)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

- n_sample:

  Number of non-zero values to sample for the preview. NULL reads all.

## Value

A list with:

- stats - named vector: min, max, mean, median, and fraction of values
  that are whole numbers

- is_integer_valued - logical; TRUE if \>99% of sampled non-zero values
  are whole numbers

- type - "CSR" or "CSC"

- dims - named integer vector c(obs, var)

- sample - numeric vector of sampled non-zero values

## Examples

``` r
# check whether /X holds raw counts before loading it
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
summary_x <- read_h5ad_x_summary(f_path)
summary_x$stats
#>               min               max              mean            median 
#>            1.0000          254.0000           12.1641            4.0000 
#> whole_number_frac 
#>            1.0000 
summary_x$is_integer_valued
#> [1] TRUE

unlink(f_path)
```
