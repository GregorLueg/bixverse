# Seurat rPCA batch correction in Rust

**\[experimental\]** This function implements the reciprocal PCA (rPCA)
anchor integration from Stuart, et al. Each batch keeps its own PCA
basis and the other batch's expression is projected into it. Anchors
live in these projected spaces and are then used to apply a
kernel-weighted correction on the union PCA embedding. Cheaper than CCA
and less aggressive in its correction.

## Usage

``` r
rs_seurat_rpca(
  f_path_gene,
  f_path_cell,
  cell_indices,
  gene_indices,
  batch_indices,
  precomputed_pca,
  rpca_params,
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

  Optional PCA matrix. If you want to provide a pre-computed matrix. The
  per-batch PCAs are always recomputed.

- rpca_params:

  List. Contains all of the Seurat rPCA parameters.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

The batch-corrected embedding space.

## References

Stuart, et al., Cell, 2019
