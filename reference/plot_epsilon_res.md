# Plot the epsilon vs. power law goodness of fit result

Plots the epsilon results (if they can be found in the class). The
x-axis reflects the different epsilon parameters for the radial basis
function, and the y-axis the R2 value that the resulting networks
follows a power law distribution (i.e., scale free topology).

## Usage

``` r
plot_epsilon_res(object)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

If epsilon results were found, returns the ggplot. Otherwise, throws a
warning and returns NULL.
