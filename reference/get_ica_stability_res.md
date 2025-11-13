# Get the ICA component data (stability, convergence, nMI)

Getter function to extract the ICA component data in terms of stability,
convergence and normalised mutual information between the components. If
not found will return `NULL`.

## Usage

``` r
get_ica_stability_res(object)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

data.table with the ICA parameter data (if found. Otherwise `NULL`.)
