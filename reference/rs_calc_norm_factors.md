# Calculate normalisation factors

**\[experimental\]** edgeR's `calcNormFactors`, via the `edge-rs` crate.

## Usage

``` r
rs_calc_norm_factors(counts, lib_size, norm_method)
```

## Arguments

- counts:

  Integer or double matrix. Raw counts of genes x samples.

- lib_size:

  Numeric vector or NULL. Library size per sample. Pass the pre-filter
  column sums after filtering genes, as edgeR keeps them. NULL uses the
  column sums of `counts`.

- norm_method:

  String. One of `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`,
  case-insensitive.

## Value

Numeric vector of normalisation factors, one per sample.

## References

Robinson and Oshlack, Genome Biol, 2010
