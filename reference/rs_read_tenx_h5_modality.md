# Loads in a modality from a 10x h5 file

**\[experimental\]** Reads one feature type (e.g. `"Antibody Capture"`)
from a 10x h5 file as a dense matrix. Needs a v3 file; v2 files carry no
feature types and error.

## Usage

``` r
rs_read_tenx_h5_modality(f_path, version, feature_type)
```

## Arguments

- f_path:

  String. The path to the h5 file.

- version:

  String. The 10x version. If `"auto"`, the version is detected from the
  file.

- feature_type:

  String. The feature type to return. Matched against the
  whitespace-trimmed feature types of the file.

## Value

A list with:

- counts - Dense numerical matrix of cells x features.

- barcodes - Character vector of the cell barcodes, in file order.

- features - Character vector of the feature names.
