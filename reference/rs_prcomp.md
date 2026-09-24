# Rust implementation of prcomp

**\[experimental\]** Runs the singular value decomposition over the
matrix x. Assumes that samples = rows, and columns = features.

## Usage

``` r
rs_prcomp(x, scale, top_pcs)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features.

- scale:

  Boolean. Shall the columns be variance normalised. (Mean centring will
  automatically occur.)

- top_pcs:

  Optional integer. Only return the top PCs (under the hood all of them
  will be calculated). `NULL` returns all.

## Value

A list with:

- scores - The product of x (centred and potentially scaled) with v.

- v - v matrix of the SVD.

- s - Standard deviations of the PCs, i.e. singular values divided by
  `sqrt(nrow(x) - 1)`.

- scaled - Boolean. Was the matrix scaled.
