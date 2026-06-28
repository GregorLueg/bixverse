# bixverse SymphonyReference class

Holds a Symphony reference: PCA loadings, per-HVG scaling stats, soft
cluster centroids and the cached compression terms (Nr, C) needed to map
queries without the reference cells. The Harmony-corrected reference
embedding (`z_corr`) is retained for downstream label transfer.
Reference cell labels (e.g. cell type annotations) can be snapshotted at
build time via `label_columns` in
[`build_symphony_ref()`](https://gregorlueg.github.io/bixverse/reference/build_symphony_ref.md)
or attached post-hoc via
[`add_symphony_labels()`](https://gregorlueg.github.io/bixverse/reference/add_symphony_labels.md).
For details on the method, refer to Kang et al.

## Usage

``` r
SymphonyReference(
  hvg_gene_names,
  gene_means,
  gene_sds,
  loadings,
  z_corr,
  centroids,
  nr,
  c_cache,
  no_pcs,
  harmony_backend,
  batch_vars,
  slim = FALSE,
  z_orig = NULL,
  r = NULL,
  labels = NULL
)
```

## Arguments

- hvg_gene_names:

  Character vector of HVG gene names in reference loading order.

- gene_means:

  Per-HVG mean of the normalised reference data.

- gene_sds:

  Per-HVG standard deviation of the normalised reference data.

- loadings:

  PCA gene loadings matrix (n_hvgs x d).

- z_corr:

  Post-Harmony corrected embedding (N x d).

- centroids:

  Cosine-normalised reference centroids (K x d).

- nr:

  Reference cluster sizes; row-sums of `r` (length K).

- c_cache:

  Cached `R * Z_corr` compression term (K x d).

- no_pcs:

  Number of principal components.

- harmony_backend:

  Which Harmony variant was used (`"v1"` or `"v2"`).

- batch_vars:

  Names of the batch variables used during reference construction.

- slim:

  Logical; if `TRUE`, `z_orig` and `r` are dropped. Default `FALSE`.

- z_orig:

  Pre-Harmony PCA scores (N x d). `NULL` in slim references.

- r:

  Soft cluster assignments (K x N). `NULL` in slim references.

- labels:

  Optional `data.table` of reference cell labels with `nrow(z_corr)`
  rows, one column per label. `NULL` if no labels stored.

## Value

Returns the `SymphonyReference` class for further operations.

## References

Kang et al., Nat. Commun., 2021
