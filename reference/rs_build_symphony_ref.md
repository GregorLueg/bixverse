# Build a Symphony reference (Rust)

**\[experimental\]** Builds the Symphony reference in Rust, see Kang et
al. Runs PCA on the HVGs, corrects it with Harmony and stores the terms
needed to map queries.

## Usage

``` r
rs_build_symphony_ref(
  f_path_gene,
  f_path_cell,
  cell_indices,
  hvg_indices,
  batch_labels,
  pca_params,
  no_pcs,
  harmony_params,
  harmony_version,
  seed,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the gene-based binary file.

- f_path_cell:

  String. Path to the cell-based binary file. Only read for the PFlogPF
  offsets if `pca_params` requests them.

- cell_indices:

  Integer vector. 0-based cell indices.

- hvg_indices:

  Integer vector. 0-based HVG indices.

- batch_labels:

  List of 0-indexed integer vectors (one per batch variable), each of
  length `cell_indices`.

- pca_params:

  List. Output of
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- no_pcs:

  Integer. Number of principal components.

- harmony_params:

  List. Output of
  [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md)
  or
  [`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md),
  matching `harmony_version`.

- harmony_version:

  String. `"v1"` or `"v2"`; anything else errors.

- seed:

  Integer. Seed for reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- gene_means - Numerical vector. Per-HVG mean of the normalised data.

- gene_sds - Numerical vector. Per-HVG standard deviation.

- loadings - Numerical matrix. PCA loadings (n_hvgs x d).

- z_orig - Numerical matrix. Pre-Harmony PCA scores (N x d).

- z_corr - Numerical matrix. Harmony-corrected embedding (N x d).

- r - Numerical matrix. Soft cluster assignments (K x N).

- centroids - Numerical matrix. Cosine-normalised centroids (K x d).

- nr - Numerical vector. Cluster sizes, the row sums of `r`.

- c - Numerical matrix. Compression term, `r` times `z_corr` (K x d).

## References

Kang et al., Nat Comm, 2021.
