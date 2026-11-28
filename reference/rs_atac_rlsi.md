# rLSI integration of an embedding over a peak or tile store

**\[experimental\]** Signac's `IntegrateEmbeddings()` with
`reduction = "rlsi"`: per-batch LSIs give the anchor space, the anchors
then correct `embedding`.

## Usage

``` r
rs_atac_rlsi(
  f_path_cells,
  f_path_genes,
  cells,
  features,
  batch,
  embedding,
  dims,
  n_comp,
  solver,
  k_anchor,
  k_score,
  k_weight,
  sd,
  l2_norm,
  knn_params,
  seed,
  verbose
)
```

## Arguments

- f_path_cells:

  String. Cell-based file.

- f_path_genes:

  String. Gene-based file.

- cells:

  Integer vector. Cells (0-based), rows of `embedding`.

- features:

  Integer vector. Features of the per-batch LSIs (0-based).

- batch:

  Integer vector. 0-based batch per cell.

- embedding:

  Numeric matrix. Cells x dimensions to integrate.

- dims:

  Integer vector. 0-based LSI components the anchors use.

- n_comp:

  Integer. Components of the per-batch LSIs.

- solver:

  String. SVD solver.

- k_anchor:

  Integer. Neighbours for the anchors.

- k_score:

  Integer. Neighbours for the anchor scores.

- k_weight:

  Integer. Anchors per cell in the correction weights.

- sd:

  Numeric. Bandwidth divisor of the weights.

- l2_norm:

  Boolean. Scale and L2 normalise the anchor spaces.

- knn_params:

  List. The kNN parameters, see
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md).

- seed:

  Integer. Seed.

- verbose:

  Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.

## Value

The integrated embedding, cells x dimensions.

## References

Stuart, et al., Nat Methods, 2021
