# SCENIC: Generating gene-regulatory networks

**\[experimental\]**

## Usage

``` r
rs_scenic_grn(
  f_path_genes,
  cell_indices,
  gene_indices,
  tf_indices,
  scenic_params,
  seed,
  verbose
)
```

## Arguments

- f_path_genes:

  Path to the `counts_genes.bin` file.

- cell_indices:

  Integer vector. 0-indexed(!) positions of cells to include in the
  analysis

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to include.

- tf_indices:

  Integer vector. 0-indexed(!) positions of the TF predictor variables
  to use in the generation of the regression learners.

- scenic_params:

  Named list. Contains all of the parameters need for SCENIC.

- seed:

  Integer. Controls reproducibility of the function.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A gene x TF importance matrix
