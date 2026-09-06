# Bixverse implementation of ssGSEA

Implementation of the bixverse version of the single sample gene set
enrichment analysis (ssGSEA), see Barbie et al.

## Usage

``` r
calc_ssgsea(exp, pathways, ssgsea_params = params_ssgsea(), .verbose = FALSE)
```

## Arguments

- exp:

  Numerical matrix. Rows represents the features, columns the
  features/genes.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- ssgsea_params:

  List. The GSVA parameters, see
  [`params_ssgsea()`](https://gregorlueg.github.io/bixverse/reference/params_ssgsea.md)
  wrapper function. This function generates a list containing:

  - alpha - Float. The exponent defining the weight of the tail in the
    random walk performed by ssGSEA.

  - min_size - Integer. Minimum size for the gene sets.

  - max_size - Integer. Maximum size for the gene sets.

  - normalise - Boolean. Shall the scores be normalised.

- .verbose:

  Boolean. Controls verbosity.

## Value

A matrix of shape pathways (that passed the thresholds) x samples.

## References

Barbie et al., Nature, 2009

## Examples

``` r
# per-sample ssGSEA scores for two gene sets
set.seed(123L)
exp_mat <- matrix(
  rnorm(200 * 10),
  nrow = 200,
  dimnames = list(sprintf("gene_%03i", 1:200), sprintf("sample_%i", 1:10))
)
pathways <- list(
  set_a = sprintf("gene_%03i", 1:20),
  set_b = sprintf("gene_%03i", 50:80)
)
round(calc_ssgsea(exp_mat, pathways)[, 1:3], 3)
#>       sample_1 sample_2 sample_3
#> set_a    0.584    0.170   -0.147
#> set_b    0.353    0.483    0.691
```
