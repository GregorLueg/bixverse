# Benjamini-Hochberg adjustment per observation across sources

decoupler's scheme: every row (observation) is adjusted on its own
across the sources.

## Usage

``` r
lm_padj_per_obs(p_values)
```

## Arguments

- p_values:

  Numerical matrix. Observations x sources.

## Value

Numerical matrix of the same shape with the adjusted p-values.
