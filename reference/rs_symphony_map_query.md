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
  index, or `NA_integer_` if absent.

- batch_labels_query:

  List of 0-indexed integer vectors (empty = no batch correction).

- params_symphony:

  Named list. Contains the parameters for the referemce generation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with z_pca, z_corr, r.

## References

Kang et al., Nat Comm, 2021.
