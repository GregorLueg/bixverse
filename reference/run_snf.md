# Run the SNF algorithm

This function will run the SNF algorithm on top of the adjacency
matrices found in the object. You can also optionally specify which
adjacency matrices to use via the `to_include` parameter.

## Usage

``` r
run_snf(object, to_include = NULL, params = NULL)
```

## Arguments

- object:

  The underlying class, see [`snf()`](snf.md).

- to_include:

  Optional string, if you wish to only use a subset of the generated
  adjacency matrices. If `NULL` all matrices will be used for the fusion
  process.

- params:

  Optional List. If you wish to overwite the already set up parameters
  for SNF, see [`params_snf()`](params_snf.md). If `NULL`, the settings
  from within the object will be used. If not NULL, the new parameters
  will be used for this modality specifically and only for this
  modality!

## Value

The class with added adjacency matrix based on the SNF algorithm.
