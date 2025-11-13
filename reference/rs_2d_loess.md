# Rust implementation of a Loess function

Rust implementation of a Loess function

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

  Numeric. The span parameter. Needs to be between 0.1 and 1.

- degree:

  Integer. Either 1 (linear) or 2 (quadratic). Other values will cause
  an error.

## Value

A list with the following items

- predicted - The predicted values.

- residuals - The residuals for every data point.

- valid_idx - Which data indices were included.
