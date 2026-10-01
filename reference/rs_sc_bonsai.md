# Build a Bonsai tree from single cell counts

**\[experimental\]** Reads the raw counts of the given genes and cells
from the binary stores, runs Sanity on them for posterior log fold
changes with error bars, builds a Bonsai tree over the cells and lays it
out in 2D. Sanity runs on the CPU.

## Usage

``` r
rs_sc_bonsai(
  f_path_gene,
  f_path_cell,
  cell_indices,
  gene_indices,
  bonsai_params,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- f_path_cell:

  String. Path to the `counts_cells.bin` file. Supplies the library
  sizes.

- cell_indices:

  Integer. The cell indices to use. (0-indexed!) Sets the leaf order.

- gene_indices:

  Integer. The gene indices to use. (0-indexed!)

- bonsai_params:

  List. Parameter list, see
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- parent - Integer. Parent of each node (0-indexed!), `-1` for the root.
  Nodes `0` to `n_leaves - 1` are the cells in `cell_indices` order.

- branch - Numeric. Branch length above each node.

- x - Numeric. Horizontal coordinate of each node.

- y - Numeric. Vertical coordinate of each node.

- n_leaves - Integer. Number of leaves.

- loglik - Numeric. Final tree loglikelihood.

- steps - List with `step` and `loglik` after each search step.

- genes_used - Integer. Genes the tree was built on. (0-indexed!)

## References

de Groot, et al., Nat Biotechnol, 2026; Breda, et al., Nat Biotechnol,
2021.
