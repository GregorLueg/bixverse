# Assemble the expression matrix a gene trend fit runs over

Returns a cells x genes matrix in `obs_cells` order. The imputed layer
is subset by name in both directions rather than positionally, since it
carries its own dimnames and a stale one would otherwise slip through
silently.

## Usage

``` r
.gene_trends_expression(object, features, use_magic, obs_cells)
```

## Arguments

- object:

  One of the single cell classes.

- features:

  Optional character vector. The genes to pull.

- use_magic:

  Boolean. Read the imputed layer rather than the counts.

- obs_cells:

  Character vector. The kept cell names, in obs order.

## Value

A numeric matrix of cells x genes.
