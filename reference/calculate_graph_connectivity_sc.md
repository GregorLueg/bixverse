# Calculate the graph connectivity per cell type

For each cell type, restricts the kNN graph to that cell type and takes
the fraction of its cells in the largest connected component. A cell
type split across batches after correction falls apart into several
components and scores low. 1 means every cell type is one connected
piece.

## Usage

``` r
calculate_graph_connectivity_sc(object, cell_type_column, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cell_type_column:

  String. The column with the cell type labels in the obs data of the
  class.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `GraphConnectivityScores` object with the following elements

- per_cell_type - Named numeric. Connectivity per cell type.

- mean_connectivity - Mean connectivity across cell types.

- median_connectivity - Median connectivity across cell types.

## References

Luecken, et al., Nat. Methods, 2022

## Examples

``` r
# connectivity of each cell type in the kNN graph
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_graph_connectivity_sc(sc, cell_type_column = "cell_grp")
#> Graph Connectivity
#>   Cell types: 3
#>   Mean:    1.0000 (1 = every cell type connected)
#>   Median:  1.0000
#>   Lowest:  cell_type_1 (1.000), cell_type_2 (1.000), cell_type_3 (1.000)

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
