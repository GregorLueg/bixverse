# Validate a k range for a consensus NMF sweep

Checked on the R side because the disk-backed single cell sweep loads
the whole count matrix before Rust gets to validate anything, so a typo
would otherwise cost a full load.

## Usage

``` r
.assert_nmf_k_range(k_range)
```

## Arguments

- k_range:

  Integer vector. The ranks to evaluate.

## Value

The range as an integer vector.
