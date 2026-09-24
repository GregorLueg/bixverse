# Calculate LISI scores on any label

**\[experimental\]** Computes the Local Inverse Simpson's Index on the
kNN graph: the effective number of labels in each cell's neighbourhood.
On batch labels this is iLISI (higher is better mixing), on cell type
labels cLISI (lower is better separation). Both come back rescaled to
`[0, 1]`, higher is better, as in scIB.

## Usage

``` r
rs_lisi(knn_mat, knn_dist, labels, perplexity, verbose)
```

## Arguments

- knn_mat:

  Integer matrix. The rows represent the cells and the columns the
  neighbour indices (0-indexed!).

- knn_dist:

  Numeric matrix or NULL. The kNN distances, same shape as `knn_mat`. If
  provided, neighbours are weighted with a perplexity-calibrated
  Gaussian kernel as in Korsunsky et al.; if NULL, neighbours are
  weighted uniformly.

- labels:

  Integer vector. The label (batch or cell type) per cell. The codes
  need not be 0-based or contiguous.

- perplexity:

  Numeric or NULL. Perplexity for the weighted version. NULL defaults to
  30; values above k are clamped to k. Ignored if `knn_dist` is NULL.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following items

- per_cell - Per-cell LISI scores

- mean_lisi - Mean LISI

- median_lisi - Median LISI

- n_labels - Number of distinct labels

- ilisi_norm - Median LISI rescaled as iLISI, `(median - 1) / (n - 1)`

- clisi_norm - Median LISI rescaled as cLISI, `(n - median) / (n - 1)`

## References

Korsunsky, et al., Nat Methods, 2019; Luecken, et al., Nat Methods, 2022
