# Seurat CCA batch correction in Rust

**\[experimental\]** This function implements the canonical correlation
analysis (CCA) anchor integration from Stuart, et al. Anchors are
identified in a per-pair L2-normalised canonical correlation embedding,
scored via shared neighbours, filtered in gene space and then used to
apply a kernel-weighted correction on the union PCA embedding.

## Usage

``` r
rs_seurat_cca(
  f_path_gene,
  f_path_cell,
  cell_indices,
  gene_indices,
  batch_indices,
  precomputed_pca,
  cca_params,
  verbose,
  seed
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- f_path_cell:

  String. Path to the `counts_cells.bin` file. Used if you wish to use
  the PFlogPF transformation during the optional PCA step.

- cell_indices:

  Integer. The cell indices to use. (0-indexed!)

- gene_indices:

  Integer. The gene indices to use. (0-indexed!) Ideally these are
  batch-aware highly variable genes.

- batch_indices:

  Integer vector. These represent to which batch a given cell belongs.
  Need to be 0-indexed and contiguous.

- precomputed_pca:

  Optional PCA matrix. If you want to provide a pre-computed matrix.

- cca_params:

  List. Contains all of the Seurat CCA parameters.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

The batch-corrected embedding space.

## References

Stuart, et al., Cell, 2019
