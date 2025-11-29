# Generate a new miloR result

This is a helper class that wraps the miloR results together. The
general idea of the approach is to use the kNN graph generated in the
single cell data, generate representative neighbourhoods and calculate
differential abundances within these neighbourhoods. For further details
on the method, please refer to Dann, et al.

## Usage

``` r
new_sc_miloR_res(nhoods, sample_counts, spatial_dist, params)
```

## Arguments

- nhoods:

  Sparse dgCMatrix with cells x neighbourhoods.

- sample_counts:

  Integer matrix. Represents neighbourhoods x cells from sample of
  interest.

- spatial_dist:

  Numeric. The spatial distance between the neighbourhoods to calculate
  the spatial FDR.

- params:

  Named list. The parameters that were used to generate these results.

## Value

An `sc_miloR` class that contains the provided data and has subsequent
methods to calculate differential abundance statistics.

## References

Dann, et al., Nat Biotechnol, 2022
