# Read in 10x h5 ADT data

Helper function to load in ADT counts from h5. Leverages Rust under the
hood for faster filtering and reading.

## Usage

``` r
read_tenx_h5_adt(f_path, feature_type = "Antibody Capture")
```

## Arguments

- f_path:

  String. File to the h5 file from which to read the ADT counts.

- feature_type:

  String. The feature type to return. Defaults here to
  `"Antibody Capture"`

## Value

A dense matrix of cells x features
