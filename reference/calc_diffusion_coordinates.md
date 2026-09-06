# Calculate diffusion coordinates

To leverage the quality metrics from Persad, et al., we need the
diffusion coordinates to then calculate if a cell is a dense or sparse
region of the manifold, its compactness and separation to other meta
cells. To do so, generate a diffusion map on the original data based on
the approach of SEACells and add the data to the object, see Persad, et
al.

## Usage

``` r
calc_diffusion_coordinates(
  object,
  knn_data,
  n_dcs = 10L,
  k_density = 150L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `MetaCells` class.

- knn_data:

  `SingleCellNearestNeighbour` class. Contains the kNN graph from the
  original cells.

- n_dcs:

  Integer. Number of diffusion coordinates to use. Defaults to `10L`.

- k_density:

  Integer. The k-th neighbour to use for the density region estimation.
  Defaults to `150L`.

- seed:

  Integer. Seed for reproducibility

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The class with the diffusion map coordinates, density distance and
region attached.

## References

Persad, et al. Nat Biotechnol, 2023

## Examples

``` r
# diffusion map off the source kNN graph, giving each meta cell a region
sc <- demo_single_cells()
mc <- generate_bt_meta_cells_sc(
  sc,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  .verbose = FALSE
)
mc <- calc_diffusion_coordinates(
  mc,
  knn_data = get_knn_obj(sc),
  .verbose = FALSE
)
table(mc[["density_region"]]$density_region)
#> 
#> high  low  mid 
#>   12   11   27 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
