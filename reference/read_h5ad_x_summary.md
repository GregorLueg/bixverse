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
