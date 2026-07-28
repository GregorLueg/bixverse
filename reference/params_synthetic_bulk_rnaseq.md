# Wrapper function to generate synthetic bulk RNAseq parameters

Parameters for
[`synthetic_bulk_cor_matrix()`](https://gregorlueg.github.io/bixverse/reference/synthetic_bulk_cor_matrix.md).
Counts come from a negative binomial with a mean-dispersion trend;
co-expression modules are planted by putting each module's genes on a
shared latent factor. The `generator` picks how loadings and factors are
drawn, which is what makes a given dataset a fair or unfair benchmark
for a given method:

- `"hub_modular"` - LogNormal loadings on a Normal factor. Some genes
  end up far more connected than others, so this is the WGCNA-style
  default.

- `"modular"` - Beta(5, 2) loadings on a Normal factor. Homogeneous
  within-module correlation and no hubs.

- `"non_negative_factor"` - LogNormal loadings on a Gamma factor. The
  activity matrix is non-negative by construction, so NMF has a ground
  truth it can actually reach.

- `"non_gaussian_factor"` - LogNormal loadings on a Laplace factor.
  Non-Gaussian sources satisfy ICA identifiability.

## Usage

``` r
params_synthetic_bulk_rnaseq(
  num_samples = 100L,
  num_genes = 1000L,
  module_sizes = c(100L, 100L, 100L),
  generator = c("hub_modular", "modular", "non_negative_factor", "non_gaussian_factor"),
  seed = 123L,
  mean_exp_gamma_shape = 5,
  mean_exp_gamma_scale = 10,
  disp_intercept = 0.2,
  disp_slope = 0.3,
  noise_std = 0.1,
  factor_std = 0.5,
  factor_shape = 2,
  factor_scale = 0.3,
  loading_mu = 0,
  loading_sigma = 0.7,
  hub_percentile = 0.1
)
```

## Arguments

- num_samples:

  Integer. Number of samples (columns) to simulate.

- num_genes:

  Integer. Number of genes (rows) to simulate.

- module_sizes:

  Integer vector. Sizes of the co-expression modules. The sum must be
  smaller or equal to `num_genes`. Genes are assigned in contiguous
  blocks from the first gene onwards, any remainder is background. Use
  `integer(0)` for no modules. Must be an integer vector, see the note
  below.

- generator:

  String. Which topology and distribution family to plant. One of
  `c("hub_modular", "modular", "non_negative_factor", "non_gaussian_factor")`,
  see the description. Defaults to `"hub_modular"`.

- seed:

  Integer. Seed for reproducibility purposes.

- mean_exp_gamma_shape, mean_exp_gamma_scale:

  Float. Shape and scale of the Gamma the per-gene mean expression is
  drawn from.

- disp_intercept, disp_slope:

  Float. Intercept and slope of the negative binomial dispersion trend
  `disp = 1 / (a + b * mean)`. This is what gives you
  heteroskedasticity: lowly expressed genes show higher variance.

- noise_std:

  Float. Per-gene per-sample noise standard deviation on the latent
  log-signal. Smaller values track the module factor more tightly and
  give stronger within-module correlation. Defaults to `0.1` rather than
  the crate's `0.3`, see the note below.

- factor_std:

  Float. Standard deviation of the Normal factor. Only used by
  `"hub_modular"` and `"modular"`; the other two generators draw their
  factor from `factor_shape`/`factor_scale` instead. Defaults to `0.5`
  rather than the crate's `0.3`, see the note below.

- factor_shape, factor_scale:

  Float. Shape and scale of the Gamma factor for
  `"non_negative_factor"`. `factor_scale` doubles as the Laplace scale
  for `"non_gaussian_factor"`.

- loading_mu, loading_sigma:

  Float. Location and scale of the LogNormal the loadings are drawn
  from. Unused by `"modular"`, which draws Beta(5, 2).

- hub_percentile:

  Float. Top fraction of module genes flagged as hubs by loading rank.
  Must be in `(0, 1]`.

## Value

A list with the parameters for usage in subsequent functions.

## Details

`noise_std` and `factor_std` default to `0.1` and `0.5`, not to the
`bixverse-rs` values of `0.3` and `0.3`. At the crate defaults the
`"modular"` generator plants modules too weakly to detect at 1000 genes
by 100 samples: the within-module minus cross-module mean absolute
Spearman gap comes out around `0.06`, against `0.17` to `0.23` for the
other three. The values here put all four generators in the `0.30` to
`0.39` band, so a comparison across generators reflects the method
rather than the signal strength it happened to be handed. Pass the crate
values explicitly if you want a harder problem.

## References

Zhang & Horvath, Stat Appl Genet Mol Biol, 2005
