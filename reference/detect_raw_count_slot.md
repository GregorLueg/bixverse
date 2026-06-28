# Detect which slot holds raw integer counts in an h5ad file

Samples the non-zero values of each candidate slot and returns the first
one that is integer-valued. Slots are checked in the given order, so
dedicated count slots take precedence over `X`.

## Usage

``` r
detect_raw_count_slot(
  f_path,
  candidates = c("layers.counts", "raw.X", "X"),
  n_sample = 10000L,
  threshold = 0.99
)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

- candidates:

  Character vector of slots to test, any of "layers.counts", "raw.X",
  "X". Order defines priority.

- n_sample:

  Number of values to sample per slot.

- threshold:

  Minimum fraction of non-zero values that must be whole numbers for a
  slot to count as raw.

## Value

The detected slot name, or NULL if none qualifies.
