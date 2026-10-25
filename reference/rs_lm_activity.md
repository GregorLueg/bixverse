# ULM / MLM activity scores (decoupler-style)

**\[experimental\]** Fits a univariate (ULM) or multivariate (MLM)
linear model per observation of the feature weights of each set against
the observed values, and returns the t-value or the coefficient as the
activity score. Observations are the columns of `mat`.

## Usage

``` r
rs_lm_activity(mat, gs_list, lm_params)
```

## Arguments

- mat:

  Numerical matrix. Features x observations. Must not contain non-finite
  values.

- gs_list:

  List. One element per set, each a list with `idx` (integer, 0-indexed
  feature positions) and optionally `weight` (numeric, aligned with
  `idx`; absent means weight 1).

- lm_params:

  List. With `method` (`"ulm"` or `"mlm"`), `stat` (`"t"` or `"coef"`)
  and `p_values` (boolean).

## Value

A list with the following elements

- es - Numerical matrix of observations x sets with the scores.

- p_values - Numerical matrix of observations x sets with the two-sided
  p-values, or `NULL` if not requested.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022.
