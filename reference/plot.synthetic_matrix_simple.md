# Plot the (simple) synthetic gene expression

Plot the (simple) synthetic gene expression

## Usage

``` r
# S3 method for class 'synthetic_matrix_simple'
plot(x, ...)
```

## Arguments

- x:

  `synthetic_matrix_simple` class. Output from
  [`synthetic_signal_matrix()`](synthetic_signal_matrix.md).

- ...:

  Additional params

## Value

A plotted heatmap showing the DEG.

## Examples

``` r
if (FALSE) { # \dontrun{

synthetic_gex <- synthetic_signal_matrix()

plot(synthetic_gex)
} # }
```
