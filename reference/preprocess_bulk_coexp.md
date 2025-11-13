# Process the raw data

Function to do general pre-processing on top of the
[`bulk_coexp()`](bulk_coexp.md). Options to do scaling, HVG selection,
etc.

## Usage

``` r
preprocess_bulk_coexp(
  object,
  hvg = NULL,
  mad_threshold = NULL,
  scaling = FALSE,
  scaling_type = c("normal", "robust"),
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see [`bulk_coexp()`](bulk_coexp.md).

- hvg:

  Integer or float. If an integer, the top `hvg` genes will be included;
  if float, the float has to be between 0 and 1, representing the
  percentage of genes to include.

- mad_threshold:

  Float. Instead of of selecting number or proportion of genes, you can
  also provide a mad_threshold.

- scaling:

  Boolean. Shall the data be scaled.

- scaling_type:

  String. You have the option to use normal scaling or robust scaling.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Returns the class with the `processed_data` data slot populated and
applied parameters added to the `params` slot.
