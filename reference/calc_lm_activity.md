# ULM / MLM activity scores (decoupler-style)

Rust implementation of decoupler's univariate (ULM) and multivariate
(MLM) linear model activity scores. Per observation, the observed values
of the features are regressed on the weights of each source's targets
(ULM: one model per source; MLM: all sources in one model), and the
t-value (or the coefficient) of the slope is the activity. Positive
means the targets move in the direction of their weights.

Behaviour worth knowing:

- Rows with non-finite values are dropped (with a warning) before
  scoring.

- Constant observations score 0 with a p-value of 1.

- MLM errors on collinear sources (duplicated or proportional target
  sets, or a set covering every feature with one weight).

- The p-values treat the features as independent.

- `fdr` is the Benjamini-Hochberg adjustment per observation across
  sources, decoupler's scheme. decoupler only adjusts ULM; here both
  methods get it, so compare `pvals` against decoupler's MLM output.

## Usage

``` r
calc_lm_activity(mat, net, lm_params = params_lm_activity())
```

## Arguments

- mat:

  Numerical matrix of features x samples with row names, or a named
  numeric vector (e.g., a contrast's t statistics), treated as one
  sample. Should be roughly normal: normalised and log-transformed
  counts, or contrast-level statistics.

- net:

  data.frame. The network in long format with the columns `source`,
  `target` and optionally `weight` (decoupler / CollecTRI layout).
  Absent weights are set to 1.

- lm_params:

  List. The parameters, see
  [`params_lm_activity()`](https://gregorlueg.github.io/bixverse/reference/params_lm_activity.md).
  This function generates a list containing:

  - method - String. `"ulm"` or `"mlm"`.

  - stat - String. `"t"` or `"coef"`.

  - p_values - Boolean. Shall p-values be computed.

  - tmin - Integer. Minimum number of targets per source in the data.

## Value

For a matrix input, a list with the following elements:

- es - Numerical matrix of samples x sources with the activities.

- pvals - Numerical matrix of samples x sources with the two-sided
  p-values, or `NULL` if `p_values = FALSE`.

- fdr - Numerical matrix of the BH-adjusted p-values (per sample across
  sources), or `NULL` if `p_values = FALSE`.

For a vector input, a data.table with the columns `source`, `es`, and,
if requested, `pvals` and `fdr`. `NULL` with a warning if no source has
at least `tmin` targets in the data.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022

## Examples

``` r
# two regulators, one with its targets shifted up in the second sample
set.seed(42L)
genes <- sprintf("gene_%03i", 1:200)
mat <- matrix(
  rnorm(400),
  ncol = 2L,
  dimnames = list(genes, c("ctrl", "treated"))
)
mat[1:15, "treated"] <- mat[1:15, "treated"] + 2
net <- data.frame(
  source = rep(c("tf_a", "tf_b"), each = 15L),
  target = genes[c(1:15, 101:115)],
  weight = 1
)
res <- calc_lm_activity(
  mat,
  net,
  lm_params = params_lm_activity(p_values = TRUE)
)
res$es
#>             tf_a        tf_b
#> ctrl    2.133174  0.03763638
#> treated 7.184559 -0.70848693
```
