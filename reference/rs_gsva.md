# Rust version of the GSVA algorithm

Rust-based implementation of the popular GSVA algorithm. Has further
performance optimisations compared to the original implementation.

## Usage

``` r
rs_gsva(exp, gs_list, tau, gaussian, max_diff, abs_rank, timings)
```

## Arguments

- exp:

  Numerical matrix. The expression matrix with rows = genes, and columns
  = samples

- gs_list:

  List. A list containing the indices of the pathway genes (needs to be
  null indexed). See [`rs_prepare_gsva_gs()`](rs_prepare_gsva_gs.md).

- tau:

  Float. Tau parameter, usual recommendation is to use `1.0` here.
  Larger values emphasise the tails more.

- gaussian:

  Boolean. If `TRUE` the Gaussian kernel will be used, if `FALSE` the
  Poisson kernel will be used.

- max_diff:

  Boolean. Scoring mode: `TRUE` = difference, `FALSE` = larger absolute
  value

- abs_rank:

  Booelan. If `TRUE` = pos-neg, `FALSE` = pos+neg

- timings:

  Boolean. Prints timings from the algorithm.

## Value

Returns a matrix of gene set ES scores x samples.
