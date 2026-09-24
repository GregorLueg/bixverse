# Restrict genes to those the fitted model covers

A gene one group filtered out never reaches the shared axis, so an HVG
set picked by another method can name genes the model has no
coefficients for. Rust would error naming a raw store index, which is
not something you can act on from R.

## Usage

``` r
.residual_gene_axis(gene_indices, fit, .verbose = TRUE)
```

## Arguments

- gene_indices:

  Integer. The 0-based genes requested.

- fit:

  `ScResidualFit` object.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The 0-based genes the model covers, ascending.
