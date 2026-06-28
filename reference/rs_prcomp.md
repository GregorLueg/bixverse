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

  Boolean. Shall the columns be variance normalised. (Mean centering
  will automatically occur.)

- top_pcs:

  Optional integer. Only return the top PCs (under the hood all of them
  will be calculated).

## Value

A list with:

- scores - The product of x (centred and potentially scaled) with v.

- v - v matrix of the SVD.

- s - Eigenvalues of the SVD.

- scaled - Boolean. Was the matrix scaled.
