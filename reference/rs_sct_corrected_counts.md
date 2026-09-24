# Writes scTransform-corrected counts to a new store

**\[experimental\]** Reverses the residual transform with every latent
variable, the library size included, held at its median, so the depth
structure is removed while the per-sample intercept is kept. The result
is written as a new gene-major store.

The output is re-indexed: its gene axis is the model's, so gene `j` in
the written store is `genes[j + 1]` of the source. It also has no
normalised layer, since corrected counts carry no library size to scale
to.

## Usage

``` r
rs_sct_corrected_counts(
  f_path_gene,
  residual_fit,
  cell_indices,
  f_path_out,
  gene_batch_size,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- residual_fit:

  List. A scTransform fit from
  [`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).
  The analytic Pearson model has no corrected-count equivalent.

- cell_indices:

  Integer vector. The cell indices to use. (0-indexed!) Must be the
  selection the fit was fitted on.

- f_path_out:

  String. Path of the gene-major file to write.

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- f_path - The file that was written.

- genes - The source gene indices of the written axis. (0-indexed!)

- n_genes - Number of genes written.

- n_cells - Number of cells written.
