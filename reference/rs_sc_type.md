# Run the ScType scoring approach

**\[experimental\]** This Rust function implements the cell type scoring
approach from Ianevski et al. (2022).

## Usage

``` r
rs_sc_type(
  f_path,
  cell_indices,
  cell_markers,
  sensitivity,
  weight_floor,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- cell_indices:

  Integer vector. 0-indexed(!) positions of cells to include in the
  analysis

- cell_markers:

  A list with the cell marker gene indices.

- sensitivity:

  Boolean. Shall a sensitivity correction be applied that downweights
  common cell type markers.

- weight_floor:

  Optional numeric. If `sensitivity = TRUE`, what is the weight floor.
  If not provided, defaults to `0.1`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with

- cell_types - String vector. The cell types

- scores - Row-major scores (cells x cell_types).

- n_cells - Number of cells

- n_cell_types - Number of cell types
