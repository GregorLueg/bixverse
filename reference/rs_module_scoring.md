# Calculate module activity scores in Rust

**\[experimental\]** Calculates module activity scores following
Seurat's `AddModuleScore`. For each module (gene set), computes the
average expression of genes in the set minus the average expression of
randomly selected control genes from the same expression bins. Genes are
binned based on their average expression across cells to ensure controls
are expression-matched.

## Usage

``` r
rs_module_scoring(
  f_path_cells,
  f_path_genes,
  gs_list,
  cells_to_keep,
  nbin,
  ctrl,
  seed,
  streaming,
  verbose
)
```

## Arguments

- f_path_cells:

  String. Path to the cell-based binary file.

- f_path_genes:

  String. Path to the gene-based binary file.

- gs_list:

  List. List of integer vectors, where each vector contains gene indices
  (0-based) for a module/gene set.

- cells_to_keep:

  Integer. Vector of indices (0-indexed) of the cells to keep.

- nbin:

  Integer. Number of bins for gene stratification.

- ctrl:

  Integer. Number of control genes to sample per gene in each module.

- seed:

  Integer. Random seed for reproducible control gene sampling.

- streaming:

  Logical. If `TRUE`, cells and genes are read in chunks to reduce
  memory usage.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

Matrix of module scores (cells x modules). Each row corresponds to a
cell from `cells_to_keep`, each column to a module from `gs_list`.

## References

Tirosh et al, Science (2016)
