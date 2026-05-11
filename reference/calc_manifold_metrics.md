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

## Value

The class with the compactness and separation scores added.

## References

Persad, et al. Nat Biotechnol, 2023
