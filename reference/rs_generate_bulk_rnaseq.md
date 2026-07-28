# Generation of bulkRNAseq-like data with optional correlation structure

**\[experimental\]** Function generates synthetic bulkRNAseq data with
heteroskedasticity (lowly expressed genes show higher variance) and
optional co-expression modules planted on a latent factor. Alongside the
counts it returns the ground truth (module membership, hub genes,
per-gene loadings and the latent factors), so downstream methods can be
scored against what was actually simulated.

## Usage

``` r
rs_generate_bulk_rnaseq(synthetic_params)
```

## Arguments

- synthetic_params:

  List. The synthetic data parameters, see
  [`params_synthetic_bulk_rnaseq()`](https://gregorlueg.github.io/bixverse/reference/params_synthetic_bulk_rnaseq.md).
  Expected elements are `num_samples`, `num_genes`, `module_sizes`
  (integer vector, empty means no modules), `generator` (one of
  `c("hub_modular", "modular", "non_negative_factor", "non_gaussian_factor")`),
  `seed`, `mean_exp_gamma_shape`, `mean_exp_gamma_scale`,
  `disp_intercept`, `disp_slope`, `noise_std`, `factor_std`,
  `factor_shape`, `factor_scale`, `loading_mu`, `loading_sigma` and
  `hub_percentile`.

## Value

List with the following elements

- counts The matrix of simulated counts. Rows are genes, columns are
  samples.

- module_membership Vector defining the module membership. `0` is
  background, `1..K` the module identifier.

- module_hubs 1-indexed positions of the genes flagged as hubs. Empty
  for the `"modular"` generator, which plants no hubs.

- loadings Per-gene loading on its module's latent factor. `0` for
  background genes.

- module_factors The latent factor matrix. Rows are modules, columns are
  samples.
