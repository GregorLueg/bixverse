# Detect which slot of an h5ad file holds raw integer counts

**\[experimental\]** Samples each candidate slot in turn and returns the
first whose non-zero values are whole numbers at a rate of at least
`threshold`. A missing slot, or one sampled to all zeros, is skipped.

## Usage

``` r
rs_h5ad_raw_slot(f_path, candidates, n_sample, threshold)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- candidates:

  Character vector. Slots to test, in priority order. Any of
  `c("layers.counts", "raw.X", "X")`.

- n_sample:

  Integer. Values to sample per slot.

- threshold:

  Numeric. Minimum fraction of whole numbers among the non-zero values.

## Value

The first qualifying entry of `candidates`, or `NULL`.
