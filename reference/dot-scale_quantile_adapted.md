# Quantile-clipped min-max scaling

Clips `x` to its `[outlier_cutoff, 1 - outlier_cutoff]` quantiles, then
maps linearly to `[0, 1]`. Returns 0.5 everywhere if the two quantiles
coincide.

## Usage

``` r
.scale_quantile_adapted(x, outlier_cutoff = 0.05)
```

## Arguments

- x:

  Numeric vector.

- outlier_cutoff:

  Quantile cutoff. Defaults to 0.05.

## Value

Numeric vector in `[0, 1]`.
