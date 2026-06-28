# Build a Symphony reference (Rust)

**\[experimental\]** Builds the Symphony reference in Rust, see Kang et
al.

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

  String. Path to the cell-based binary file.

- cell_indices:

  Integer vector. 0-based cell indices.

- hvg_indices:

  Integer vector. 0-based HVG indices.

- batch_labels:

  List of 0-indexed integer vectors (one per batch variable).

- pca_params:

  List. Output of
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- no_pcs:

  Integer.

- harmony_params:

  List. Output of
  [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md)
  or
  [`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md).

- harmony_version:

  String. "v1" or "v2".

- seed:

  Integer.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with gene_means, gene_sds, loadings, z_orig, z_corr, r,
centroids, nr, c.

## References

Kang et al., Nat Comm, 2021.
