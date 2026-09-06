# Get the offsets for the CLR/PFlogPF transformation prior PCA

Helper function to get the offsets for the CLR/PFlogPF transformation.

## Usage

``` r
mc_get_clr_offsets(object, cell_indices = NULL)
```

## Arguments

- object:

  `MetaCells` class.

- cell_indices:

  Optional integer. Defines the indices of the (meta)cells to use for
  the calculation.

## Value

A vector of length cell_indices which contains the CLR offsets

## Examples

``` r
# one offset per meta cell
sc <- demo_single_cells()
mc <- generate_bt_meta_cells_sc(
  sc,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  .verbose = FALSE
)
head(mc_get_clr_offsets(mc))
#> meta_cell_01 meta_cell_02 meta_cell_03 meta_cell_04 meta_cell_05 meta_cell_06 
#>   0.01936928   0.01946831   0.01936383   0.01938593   0.01938824   0.01945201 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
