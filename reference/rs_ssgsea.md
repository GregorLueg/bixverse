# Rust version of the ssGSEA algorithm

Rust-based implementation of the popular single sample GSEA algorithm.
Has further performance optimisations compared to the original
implementation.

## Usage

``` r
rs_ssgsea(exp, gs_list, alpha, normalise, timings)
```

## Arguments

- exp:

  Numerical matrix. The expression matrix with rows = genes, and columns
  = samples

- gs_list:

  List. A list containing the indices of the pathway genes (needs to be
  null indexed). See [`rs_prepare_gsva_gs()`](rs_prepare_gsva_gs.md).

- alpha:

  Float. The alpha parameter to adjust the weights.

- normalise:

  Boolean. Shall the scores be normalised.

- timings:

  Boolean. Prints timings from the algorithm.

## Value

Returns a matrix of gene set ES scores x samples.
