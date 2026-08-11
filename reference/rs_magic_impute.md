# Impute a subset of genes with MAGIC

**\[experimental\]** Implementation of MAGIC in Rust. A row-stochastic
diffusion operator is built over the provided kNN graph and applied
`n_steps` times to the counts of the requested genes. The operator
matches the one Palantir builds, i.e. an adaptive bandwidth taken from
the `k/3`-th neighbour distance, which differs from the graphtools
kernel the reference implementation uses.

Deliberately restricted to a gene subset. The output is dense, and
smoothing over overlapping neighbourhoods inflates gene-gene
correlation, so imputed counts must not be fed into correlation-based
methods such as Hotspot, SCENIC, differential correlation or CoReMo.

## Usage

``` r
rs_magic_impute(
  f_path,
  knn_data,
  cell_indices,
  total_cells,
  gene_indices,
  magic_params,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- knn_data:

  List. The `SingleCellNearestNeighbour` data with `indices`
  (0-indexed!), `dist`, `k` and `dist_metric`. The indices are positions
  within `cell_indices`, not global cell ids. Whether the distances are
  treated as squared is derived from `dist_metric`.

- cell_indices:

  Integer vector. The global cell indices (0-indexed!) the kNN graph was
  built over, in kNN row order.

- total_cells:

  Integer. The cell count of the binary store, not of the selection.

- gene_indices:

  Integer vector. Gene indices (0-indexed!) to impute.

- magic_params:

  List. Parameter list, see
  [`params_sc_magic()`](https://gregorlueg.github.io/bixverse/reference/params_sc_magic.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

Numerical matrix of cells x genes with the imputed counts. Rows follow
`cell_indices`, columns follow `gene_indices`.

## References

van Dijk, et al., Cell, 2018.
