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

  Float. Tau parameter, usual recommendation is to use `1.0` here.
  Larger values emphasise the tails more.

- min_size:

  Integer. Minimum number of genes per gene set.

- max_size:

  Integer. Maximum number of genes per gene set.

- max_diff:

  Boolean. Scoring mode for GSVA, if `TRUE` = difference; if `FALSE` =
  larger absolute value.

- abs_rank:

  Boolean. If `TRUE` = pos - neg, `FALSE` = pos + neg.

## Value

List with parameters for usage in subsequent function.
