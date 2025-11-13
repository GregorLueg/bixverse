# Helper function to calculate the induced sparsity

Helper function to calculate the induced sparsity

## Usage

``` r
calculate_sparsity_stats(object, no_exp_bins = 10L)
```

## Arguments

- object:

  `synthetic_bulk_data` object. You need to have run
  [`simulate_dropouts()`](simulate_dropouts.md) for this function to
  work.

- no_exp_bins:

  Integer. Number of expression bins to check. Defaults to `10L`.

## Value

A list with various statistics about the sparsity

- original_sparsity - Original proportion of zeroes in the counts.

- final_sparsity - Sparsity after applying
  [`simulate_dropouts()`](simulate_dropouts.md).

- added_sparsity - Added sparsity.

- gene_sparsity_mean - Mean sparsity for the genes.

- gene_sparsity_sd - SD sparsity form the genes.

- sample_sparsity_mean - Mean sparsity for the genes.

- sample_sparsity_sd - SD sparsity for the genes.

- dropout_by_expression - Dropout per expression bin level.
