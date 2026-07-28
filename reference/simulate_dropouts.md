# Simulate sequencing-depth dropouts on synthetic bulk data

This function induces Splatter-style sequencing-depth sparsity on the
data. Per sample a size factor
`s_j ~ LogNormal(0, capture_efficiency_sigma)` is drawn, giving a target
library size of `target_library_size * s_j`. Each gene is then
binomially thinned to approach that target, so dropout falls out of the
library size rather than an explicit per-gene dropout curve. Retention
probability is capped at 1, meaning samples already below their target
are left alone rather than upsampled.

## Usage

``` r
simulate_dropouts(object, sparsity_params = params_bulk_sparsity())
```

## Arguments

- object:

  The `synthetic_bulk_data` class.

- sparsity_params:

  List. The sparsification parameters, see
  [`params_bulk_sparsity()`](https://gregorlueg.github.io/bixverse/reference/params_bulk_sparsity.md).

## Value

`synthetic_bulk_data` with added sparse data.

## References

Zappia, et al., Genome Biol, 2017
