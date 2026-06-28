# Loads in a modality from a 10x h5 file

**\[experimental\]**

## Usage

``` r
rs_read_tenx_h5_modality(f_path, version, feature_type)
```

## Arguments

- f_path:

  String. The path to the h5 file

- version:

  String. The 10x version. If `"auto"` uses the automatic detection.

- feature_type:

  String. The feature type to return.

## Value

A list with:

- counts - Numerical matrix of cells x features

- barcodes - The barcodes as a string

- features - The features as a string
