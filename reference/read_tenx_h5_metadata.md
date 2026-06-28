# Read barcode and feature tables and metadata from a 10x h5 file

Useful for exploring the data stored in a 10x CellRanger v2/v3 h5 file,
including the breakdown of feature types in a multi-modal file.

## Usage

``` r
read_tenx_h5_metadata(f_path)
```

## Arguments

- f_path:

  File path to the 10x `.h5` file.

## Value

A list with:

- obs - data.table of barcodes

- var - data.table of features (id, name, and feature_type for v3)

- dims - named integer vector c(obs, var)

- version - "v2" or "v3"

- feature_types - named integer vector of feature_type counts (v3 only,
  otherwise NULL)
