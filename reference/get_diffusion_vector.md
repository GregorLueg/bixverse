# Get the diffusion vector

Returns the diffusion vector if you ran
[`tied_diffusion()`](https://gregorlueg.github.io/bixverse/reference/tied_diffusion.md)
or
[`diffuse_seed_nodes()`](https://gregorlueg.github.io/bixverse/reference/diffuse_seed_nodes.md).

## Usage

``` r
get_diffusion_vector(object)
```

## Arguments

- object:

  The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

## Value

The diffusion vector if found. If you did not run either diffusion
functions, it will return `NULL` and a warning.
