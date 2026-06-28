# FastMNN batch correction in Rust

**\[experimental\]** This function implements the (fast) MNN algorithm
from Haghverdi, et al. Instead of working on the full matrix, it uses
under the hood PCA and generates a batch-aligned embedding space.

## Usage

``` r
rs_mnn(
  f_path_gene,
  f_path_cell,
  cell_indices,
  gene_indices,
  batch_indices,
  precomputed_pca,
  mnn_params,
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

- precomputed_pca:

  Optional PCA matrix. If you want to provide a pre-computed matrix.

- mnn_params:

  List. Contains all of the fastMNN parameters.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

The batch-corrected embedding space.
