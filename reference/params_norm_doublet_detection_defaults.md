# Helper function to generate normalisation defaults for doublet detection.

Helper function to generate normalisation defaults for doublet
detection.

## Usage

``` r
params_norm_doublet_detection_defaults()
```

## Value

A list with the following parameters for normalisation specifically
designed for doublet detection methods:

- log_transform - Boolean. Shall the counts be log-normalised. Defaults
  to `TRUE`.

- mean_center - Boolean. Shall mean centreing be applied. Defaults to
  `FALSE`.

- normalise_variance - Boolean. Shall the variance be normalised.
  Defaults to `FALSE`.

- target_size - Target library size. Defaults to `1e6`
