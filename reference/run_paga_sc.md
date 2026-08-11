# Run PAGA graph abstraction

PAGA abstracts the single cell kNN graph into a cluster-level graph,
where the edge weight between two clusters measures how much more
connected they are than expected under a random null model. It also
returns the maximum spanning forest of that graph, which is the backbone
typically plotted. For details, please refer to Wolf, et al.

Only the kNN indices are used, no distances and no expression data, so
the method is safe on any of the available graphs, including the WNN
one. Run
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md)
and a clustering first.

## Usage

``` r
run_paga_sc(
  object,
  cluster_col,
  modality = c("rna", "adt", "wnn"),
  .verbose = TRUE
)
```

## Arguments

- object:

  One of `SingleCells`, `SingleCellsSubset`, `MetaCells` or
  `SingleCellsMultiModal`.

- cluster_col:

  String. The obs column holding the cluster assignments. Empty factor
  levels are retained.

- modality:

  String. One of `c("rna", "adt", "wnn")`. Which kNN graph to run over.
  Anything but `"rna"` requires a `SingleCellsMultiModal` object.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `PagaRes` S3 object with:

- connectivities - Symmetric sparse matrix of clusters x clusters with a
  zero diagonal. Values lie in `(0, 1]`.

- connectivities_tree - Maximum spanning forest of `connectivities`,
  carrying the original connectivity values on the retained edges.

- sizes - data.table with `cluster` and `n_cells`.

- params - List with `cluster_col` and `modality`.

## References

Wolf, et al., Genome Biol., 2019.
