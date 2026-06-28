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
