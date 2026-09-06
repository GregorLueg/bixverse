# Helper function to calculate the induced sparsity

Helper function to calculate the induced sparsity

## Usage

``` r
calculate_sparsity_stats(object, no_exp_bins = 10L)
```

## Arguments

- object:

  `synthetic_bulk_data` object. You need to have run
  [`simulate_dropouts()`](https://gregorlueg.github.io/bixverse/reference/simulate_dropouts.md)
  for this function to work.

- no_exp_bins:

  Integer. Number of expression bins to check. Defaults to `10L`.

## Value

A list with various statistics about the sparsity

- original_sparsity - Original proportion of zeroes in the counts.

- final_sparsity - Sparsity after applying
  [`simulate_dropouts()`](https://gregorlueg.github.io/bixverse/reference/simulate_dropouts.md).

- added_sparsity - Added sparsity.

- gene_sparsity_mean - Mean sparsity for the genes.

- gene_sparsity_sd - SD sparsity form the genes.

- sample_sparsity_mean - Mean sparsity for the genes.

- sample_sparsity_sd - SD sparsity for the genes.

- dropout_by_expression - Dropout per expression bin level.

## Examples

``` r
# how much sparsity the dropout simulation actually added
syn <- synthetic_bulk_cor_matrix()
syn <- simulate_dropouts(syn, params_bulk_sparsity())
stats <- calculate_sparsity_stats(syn)
stats$added_sparsity
#> [1] 0.00743
stats$dropout_by_expression
#> (0.685,1.51]  (1.51,2.33]  (2.33,3.15]  (3.15,3.97]  (3.97,4.79]  (4.79,5.61] 
#> 3.678290e-01 8.080464e-02 1.103286e-02 8.630744e-04 2.717613e-05 0.000000e+00 
#>  (5.61,6.43]  (6.43,7.25]  (7.25,8.08]   (8.08,8.9] 
#> 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
```
