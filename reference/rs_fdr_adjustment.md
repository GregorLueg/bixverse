# Calculate a BH-based FDR

Rust implementation that will be faster if you have an terrifying amount
of p-values to adjust.

## Usage

``` r
rs_fdr_adjustment(pvals)
```

## Arguments

- pvals:

  Numeric vector. The p-values you wish to adjust.

## Value

The Benjamini-Hochberg adjusted p-values.
