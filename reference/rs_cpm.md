# Counts per million

**\[experimental\]** edgeR's `cpm` on a plain count matrix, via the
`edge-rs` crate.

## Usage

``` r
rs_cpm(counts, lib_size, log, prior_count)
```

## Arguments

- counts:

  Integer or double matrix. Raw counts of genes x samples.

- lib_size:

  Numeric vector or NULL. Library size per sample, e.g.
  `lib.size * norm.factors`. NULL uses the column sums.

- log:

  Boolean. Return log2-CPM.

- prior_count:

  Numeric. Prior count added before the log. Ignored if `log = FALSE`.

## Value

Numeric matrix of (log2-)CPM values, genes x samples.
