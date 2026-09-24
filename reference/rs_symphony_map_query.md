# Map a query onto a Symphony reference (Rust)

**\[experimental\]** Maps a given other SingleCells data set to the
Symphony reference, see Kang et al.

## Usage

``` r
rs_symphony_map_query(
  f_path_query,
  cell_indices_query,
  gene_means,
  gene_sds,
  loadings,
  centroids,
  nr,
  c_cache,
  ref_to_query_gene_map,
  batch_labels_query,
  params_symphony,
  verbose
)
```

## Arguments

- f_path_query:

  String. Path to the query gene-based binary file.

- cell_indices_query:

  Integer vector. 0-based query cell indices.

- gene_means, gene_sds:

  Numerical vectors. Reference per-HVG stats.

- loadings:

  Reference PCA loadings (n_hvgs x d).

- centroids:

  Reference centroids (K x d).

- nr:

  Reference cluster sizes (length K).

- c_cache:

  Reference compression term R\*Z_corr (K x d).

- ref_to_query_gene_map:

  Integer vector. For each reference HVG slot, the 0-based query gene
  index, or `NA_integer_` (or any negative value) if absent. Absent
  slots are filled with zeros.

- batch_labels_query:

  List of 0-indexed integer vectors, one per batch variable. An empty
  list skips the batch correction (`z_corr = z_pca`).

- params_symphony:

  Named list. The query mapping parameters.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- z_pca - Numerical matrix. Query projected into the reference PC space
  (N_q x d).

- z_corr - Numerical matrix. Query after the batch correction (N_q x d).

- r - Numerical matrix. Query soft assignments onto the reference
  centroids (K x N_q).

## References

Kang et al., Nat Comm, 2021.
