# Wrapper function to generate ULM / MLM activity parameters

Wrapper function to generate ULM / MLM activity parameters

## Usage

``` r
params_lm_activity(
  method = c("ulm", "mlm"),
  stat = c("t", "coef"),
  p_values = FALSE,
  tmin = 5L
)
```

## Arguments

- method:

  String. Univariate (one model per set) or multivariate (all sets in
  one model) linear model. One of `c("ulm", "mlm")`. Defaults to
  `"ulm"`.

- stat:

  String. Return the t-value or the coefficient as the activity score.
  One of `c("t", "coef")`. Defaults to `"t"`.

- p_values:

  Boolean. Shall the two-sided p-values be computed. Defaults to
  `FALSE`.

- tmin:

  Integer. Minimum number of targets per source after removing targets
  absent from the data. Defaults to `5L`.

## Value

A named list with the following elements:

- method - String. Univariate (one model per set) or multivariate (all
  sets in one model) linear model. One of `c("ulm", "mlm")`. Defaults to
  `"ulm"`.

- stat - String. Return the t-value or the coefficient as the activity
  score. One of `c("t", "coef")`. Defaults to `"t"`.

- p_values - Boolean. Shall the two-sided p-values be computed. Defaults
  to `FALSE`.

- tmin - Integer. Minimum number of targets per source after removing
  targets absent from the data. Defaults to `5L`.

## References

Badia-i-Mompel, et al., Bioinform. Adv., 2022
