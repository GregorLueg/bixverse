# Generate a new miloR result

This is a helper class that wraps the miloR results together. The
general idea of the approach is to use the kNN graph generated in the
single cell data, generate representative neighbourhoods and calculate
differential abundances within these neighbourhoods. For further details
on the method, please refer to Dann, et al.

## Usage

``` r
new_sc_miloR_res(nhoods, sample_counts, spatial_dist, nhood_overlap, params)
```

## Arguments

- nhoods:

  Sparse dgCMatrix with cells x neighbourhoods.

- sample_counts:

  Numeric matrix. Represents neighbourhoods x samples, i.e. the cells of
  each sample found in each neighbourhood.

- spatial_dist:

  Numeric. The distance to the k-th nearest neighbour per neighbourhood,
  the `"k-distance"` weighting for the spatial FDR.

- nhood_overlap:

  Numeric. The cells each neighbourhood shares with all the others, the
  `"graph-overlap"` weighting for the spatial FDR. Taken at construction
  because it is a function of the neighbourhood matrix alone.

- params:

  Named list. The parameters that were used to generate these results.

## Value

An `miloR` class that contains the provided data and has subsequent
methods to calculate differential abundance statistics.

## References

Dann, et al., Nat Biotechnol, 2022
