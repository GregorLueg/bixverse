# Calculate manifold metrics

This function will calculate the compactness and separation of your
metacells on the manifold (defined by the diffusion map). You must have
run
[`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md)
before calling this function. The idea is that compactness indicates how
tight the metacell spans the manifold, whereas separation indicates how
well the different metacells span the manifold.

## Usage

``` r
calc_manifold_metrics(object)
```

## Arguments

- object:

  `MetaCells` class for which to calculate the different metrics.

## Value

The class with the compactness and separation scores added.

## References

Persad, et al. Nat Biotechnol, 2023

## Examples

``` r
# compactness and separation on top of the diffusion coordinates
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
mc <- calc_manifold_metrics(mc)
head(mc[[c("compactness", "separation")]])
#>    compactness separation
#>          <num>      <num>
#> 1:  0.02674035 0.17462617
#> 2:  0.03287052 0.15380448
#> 3:  0.01897743 0.22679596
#> 4:  0.04910118 0.27503037
#> 5:  0.03481250 0.04071459
#> 6:  0.01309141 0.17024054

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
