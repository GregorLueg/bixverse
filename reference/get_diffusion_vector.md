# Get the diffusion vector

Returns the diffusion vector if you ran
[`tied_diffusion()`](tied_diffusion.md) or
[`diffuse_seed_nodes()`](diffuse_seed_nodes.md).

## Usage

``` r
get_diffusion_vector(object)
```

## Arguments

- object:

  The underlying class [`network_diffusions()`](network_diffusions.md).

## Value

The diffusion vector if found. If you did not run either diffusion
functions, it will return `NULL` and a warning.
