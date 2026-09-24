# Residual variance and the variable features it selects

**\[experimental\]** Computes the per-gene residual variance from a
fitted model and selects the variable features from it. With more than
one group the selection follows Seurat: rank within each group, take the
top `n_hvg` of each, and union them, so a marker only one sample carries
is not buried by a pooled ranking. The returned set can therefore be
larger than `n_hvg`.

Each gene's residual row is regenerated, reduced and dropped, so memory
is one row per worker rather than a genes-by-cells matrix.

## Usage

``` r
rs_sc_residual_variance(
  f_path_gene,
  residual_fit,
  cell_indices,
  n_hvg,
  gene_batch_size,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- residual_fit:

  List. A fit from
  [`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).

- cell_indices:

  Integer vector. The cell indices to use. (0-indexed!) Must be the
  selection the fit was fitted on.

- n_hvg:

  Integer. Variable features to take from each group.

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- genes - The gene indices the variances are indexed by. (0-indexed!)

- variance - Matrix of residual variance, genes by groups.

- hvg - The selected gene indices, ascending. (0-indexed!)

## References

Seurat v5, `SCTransform.StdAssay`
