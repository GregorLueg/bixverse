# Wrapper function to generate GSVA parameters

Wrapper function to generate GSVA parameters

## Usage

``` r
params_gsva(
  tau = 1,
  min_size = 5L,
  max_size = 500L,
  max_diff = TRUE,
  abs_rank = FALSE
)
```

## Arguments

- tau:

  Numeric. Tau parameter, usual recommendation is to use `1.0` here.
  Larger values emphasise the tails more. Defaults to `1.0`.

- min_size:

  Integer. Minimum number of genes per gene set. Defaults to `5L`.

- max_size:

  Integer. Maximum number of genes per gene set. Defaults to `500L`.

- max_diff:

  Boolean. Scoring mode for GSVA, if `TRUE` = difference; if `FALSE` =
  larger absolute value. Defaults to `TRUE`.

- abs_rank:

  Boolean. If `TRUE` = pos - neg, `FALSE` = pos + neg. Defaults to
  `FALSE`.

## Value

A named list with the following elements:

- tau - Numeric. Tau parameter, usual recommendation is to use `1.0`
  here. Larger values emphasise the tails more. Defaults to `1.0`.

- min_size - Integer. Minimum number of genes per gene set. Defaults to
  `5L`.

- max_size - Integer. Maximum number of genes per gene set. Defaults to
  `500L`.

- max_diff - Boolean. Scoring mode for GSVA, if `TRUE` = difference; if
  `FALSE` = larger absolute value. Defaults to `TRUE`.

- abs_rank - Boolean. If `TRUE` = pos - neg, `FALSE` = pos + neg.
  Defaults to `FALSE`.
