# Voom-transform a count matrix

**\[experimental\]** limma's `voom` on counts that are already filtered:
log2-CPM against the supplied (effective) library sizes, the
mean-variance trend and the precision weights. No filtering and no
normalisation happen in here; pass `lib.size * norm.factors` as
`lib_size` to get voom on a normalised DGEList.

## Usage

``` r
rs_voom_normalise(counts, design, lib_size, span, adaptive_span)
```

## Arguments

- counts:

  Integer or double matrix. Raw counts of genes x samples.

- design:

  Numeric matrix. The design matrix of samples x coefficients. Must be
  full rank.

- lib_size:

  Numeric vector. The effective library size per sample.

- span:

  Numeric. Lowess span, only used if `adaptive_span = FALSE`.

- adaptive_span:

  Boolean. Derive the span from the number of genes, as limma does since
  3.56.

## Value

A list with the following elements

- e - Numeric matrix. The log2-CPM values, genes x samples. limma's `E`.

- weights - Numeric matrix. The precision weights, genes x samples.

- trend_x - The mean-variance trend abscissae (average log2 count).

- trend_y - The mean-variance trend ordinates (sqrt standard deviation).

- amean - Average log2-CPM per gene.

## References

Law, et al., Genome Biol, 2014
