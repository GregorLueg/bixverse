# Weighted Benjamini-Hochberg over overlapping neighbourhoods

**\[experimental\]** Milo's spatial FDR. Neighbourhoods overlap, so the
tests are not independent and a plain BH is anti-conservative. Each
p-value is weighted by the reciprocal of its connectivity and the
step-up runs on those weights. Non-finite p-values are carried through
untouched and take no part in the adjustment.

## Usage

``` r
rs_spatial_fdr(p_values, connectivity)
```

## Arguments

- p_values:

  Numeric vector. One raw p-value per tested neighbourhood.

- connectivity:

  Numeric vector. The matching connectivity per neighbourhood, either
  the k-th neighbour distances or the neighbourhood overlaps. A zero
  connectivity gets a weight of one, as in the upstream.

## Value

The adjusted p-values, in the input order.

## References

Dann, et al., Nat Biotechnol, 2022
