# Add a data modality for SNF generation

This function will add a modality as an affinity matrix for subsequent
fusion to the object. To note: all of the modalities need to have the
same number of rows, i.e., samples!

## Usage

``` r
add_snf_data_modality(object, data, data_name, params = NULL)
```

## Arguments

- object:

  The underlying class, see [`snf()`](snf.md).

- data:

  matrix or data.table. The data to transform into adjacency data and
  add to the class. Any data supplied will be assumed to be samples x
  features. The provided data type can be a data.table (for categorical
  and/or mixed types) or a matrix (for continous types). If you provide
  a data.table, the function will assume the first column are the sample
  identifiers. Please ensure that setting.

- data_name:

  String. The data modality name.

- params:

  Optional List. If you wish to overwite the already set up parameters
  for SNF, see [`params_snf()`](params_snf.md). If `NULL`, the settings
  from within the object will be used. If not NULL, the new parameters
  will be used for this modality specifically and only for this
  modality!

## Value

The class with added adjacency matrix for this data.
