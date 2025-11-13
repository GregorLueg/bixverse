# FastMNN batch correction in Rust

This function implements the (fast) MNN algorithm from Haghverdi, et al.
Instead of working on the full matrix, it uses under the hood PCA and
generates an aligned embedding space.

## Usage

``` r
rs_mnn(
  f_path_gene,
  cell_indices,
  gene_indices,
  batch_indices,
  mnn_params,
  verbose,
  seed
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- cell_indices:

  Integer. The cell indices to use. (0-indexed!)

- gene_indices:

  Integer. The gene indices to use. (0-indexed!) Ideally these are
  batch-aware highly variable genes.

- batch_indices:

  Integer vector. These represent to which batch a given cell belongs.

- mnn_params:

  List. Contains all of the fastMNN parameters.

- verbose:

  Boolean. Controls verbosity of the function.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

The batch-corrected embedding space.
