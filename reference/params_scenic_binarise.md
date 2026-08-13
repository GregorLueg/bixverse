# Wrapper function for parameters for the SCENIC binarisation

Each regulon gets its own threshold. A two-component Gaussian mixture is
fitted and compared against a single Gaussian by BIC; if the mixture
wins, the threshold is the kernel density minimum between the two
component means, otherwise it falls back to `mean + 2 * sd`. This
follows pySCENIC. AUCell fits six candidates and then lets the density
trough override all of them whenever one exists, so the two land in much
the same place.

Turn `bw_adjust` up if shallow wobbles in the density are being picked
up as troughs. AUCell effectively runs at `2`.

## Usage

``` r
params_scenic_binarise(bw_adjust = 1, n_grid = 512L, n_bins = 512L)
```

## Arguments

- bw_adjust:

  Float. Multiplier on the Silverman bandwidth of the kernel density
  estimate. Higher values smooth more. Defaults to `1`.

- n_grid:

  Integer. Number of points at which the density is evaluated between
  the two component means. Defaults to `512L`.

- n_bins:

  Integer. Number of histogram bins used to approximate the density.
  Defaults to `512L`.

## Value

A list with the binarisation parameters.

## References

Aibar, et al., Nat Methods, 2017
