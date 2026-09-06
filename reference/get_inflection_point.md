# Identify the inflection point for elbow-like data

This function will identify the index of the inflection point of an
arbitrary series `x ~ y` via the biggest increase in the first
derivative. Useful to identify key points in elbow plots.

## Usage

``` r
get_inflection_point(x, y, span = 0.5)
```

## Arguments

- x, y:

  The x and y values.

- span:

  The span parameter for the loess function.

## Value

A list containing:

- inflection_idx - Index of the inflection point

- gradient_change - Absolute change in the first derivative

## Examples

``` r
# elbow of a decaying curve that flattens out
x <- 1:20
y <- c(exp(-x[1:10] / 2), rep(0.005, 10))
get_inflection_point(x, y)$inflection_idx
#> [1] 3
```
