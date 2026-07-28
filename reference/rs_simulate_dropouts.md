# Sparsify bulkRNAseq like data

**\[experimental\]** This function takes in a (raw) count matrix (for
example from the synthetic data in bixverse) and applies Splatter-style
sequencing-depth dropout to it. Per sample a size factor
`s_j ~ LogNormal(0, capture_efficiency_sigma)` is drawn, giving a target
library size of `target_library_size * s_j`. Each gene is then
binomially thinned to approach that target. Retention probability is
capped at 1, so samples below their target are left alone rather than
upsampled.

## Usage

``` r
rs_simulate_dropouts(count_mat, sparsity_params)
```

## Arguments

- count_mat:

  Numerical matrix. Original numeric matrix. Rows are genes, columns are
  samples.

- sparsity_params:

  List. The sparsity parameters, see
  [`params_bulk_sparsity()`](https://gregorlueg.github.io/bixverse/reference/params_bulk_sparsity.md).
  Expected elements are `strategy`, `target_library_size`,
  `capture_efficiency_sigma` and `seed`.

## Value

The sparsified matrix based on the provided parameters.

## References

Zappia, et al., Genome Biol, 2017
