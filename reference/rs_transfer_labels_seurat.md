# Seurat label transfer

**\[experimental\]** Seurat's `FindTransferAnchors(reduction = "cca")`
followed by `TransferData()`. The anchors come from a CCA between the
reference and the query on the given features; the label scores are
weighted by the anchors' neighbourhoods in `weight_embedding`. Works for
an RNA or a gene activity query.

## Usage

``` r
rs_transfer_labels_seurat(
  ref_genes,
  ref_cells,
  ref_features,
  query_genes,
  query_cells,
  query_features,
  ref_labels,
  n_labels,
  weight_embedding,
  n_cc,
  k_anchor,
  k_score,
  scale_max,
  l2_norm,
  k_weight,
  sd,
  knn_params,
  seed,
  verbose
)
```

## Arguments

- ref_genes:

  String. Path to the gene-based file of the reference.

- ref_cells:

  Integer vector. Reference cells (0-indexed).

- ref_features:

  Integer vector. Reference features (0-indexed).

- query_genes:

  String. Path to the gene-based file of the query.

- query_cells:

  Integer vector. Query cells (0-indexed), rows of `weight_embedding`.

- query_features:

  Integer vector. Query features (0-indexed), matching `ref_features` in
  order.

- ref_labels:

  Integer vector. 0-indexed label code per reference cell.

- n_labels:

  Integer. Number of label codes.

- weight_embedding:

  Numeric matrix. Query cells x dimensions, the embedding of the
  weights.

- n_cc:

  Integer. Number of canonical correlation vectors.

- k_anchor:

  Integer. Neighbours for the anchors.

- k_score:

  Integer. Neighbours for the anchor scores.

- scale_max:

  Numeric. Clip of the scaled features.

- l2_norm:

  Boolean. L2 normalise the CCA space.

- k_weight:

  Integer. Anchors per query cell in the weights.

- sd:

  Numeric. Bandwidth divisor of the weights.

- knn_params:

  List. The kNN parameters, see
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md).

- seed:

  Integer. Seed.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

Numeric matrix of query cells x labels with the label scores.

## References

Stuart, et al., Cell, 2019
