# Detect which slot holds raw integer counts in an h5ad file

Samples the non-zero values of each candidate slot and returns the first
one that is integer-valued. Slots are checked in the given order, so
dedicated count slots take precedence over `X`.

## Usage

``` r
detect_raw_count_slot(
  f_path,
  candidates = c("layers.counts", "raw.X", "X"),
  n_sample = 10000L,
  threshold = 0.99
)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

- candidates:

  Character vector of slots to test, any of "layers.counts", "raw.X",
  "X". Order defines priority.

- n_sample:

  Number of values to sample per slot.

- threshold:

  Minimum fraction of non-zero values that must be whole numbers for a
  slot to count as raw.

## Value

The detected slot name, or NULL if none qualifies.

## Examples

``` r
# the synthetic writer only fills /X, and it holds raw counts
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
detect_raw_count_slot(f_path)
#> [1] "X"

unlink(f_path)
```
