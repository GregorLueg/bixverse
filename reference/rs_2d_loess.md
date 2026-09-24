# Rust implementation of a Loess function

**\[experimental\]** Fits a Loess regression of `y` on `x`. Points where
either value is non-finite are dropped from the fit.

## Usage

``` r
rs_2d_loess(x, y, span, degree)
```

## Arguments

- x:

  Numeric. The x values to fit.

- y:

  Numeric. The y values to fit.

- span:

  Numeric. The span parameter. Needs to be in `(0, 1]`.

- degree:

  Integer. Either 1 (linear) or 2 (quadratic). Other values will cause
  an error.

## Value

A list with the following items

- predicted - The predicted values, `0` for dropped points.

- residuals - The residuals for every data point, `0` for dropped
  points.

- valid_idx - 1-based indices of the points included in the fit.
