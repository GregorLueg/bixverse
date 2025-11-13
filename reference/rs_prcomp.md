# Rust implementation of prcomp

Runs the singular value decomposition over the matrix x. Assumes that
samples = rows, and columns = features.

## Usage

``` r
rs_prcomp(x, scale)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features.

- scale:

  Boolean. Shall the columns additionally be scaled.

## Value

A list with:

- scores - The product of x (centred and potentially scaled) with v.

- v - v matrix of the SVD.

- s - Eigenvalues of the SVD.

- scaled - Boolean. Was the matrix scaled.
