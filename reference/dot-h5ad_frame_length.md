# Number of rows of an h5ad dataframe group

Reads the length of the index named by the `_index` attribute. The index
is a plain dataset for older writers, or a group holding `values`
(nullable strings, as written with pandas \>= 3) or `codes`
(categoricals). Files without an `_index` attribute fall back to the
first plain dataset in the group.

## Usage

``` r
.h5ad_frame_length(f_path, h5_content, group)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

- h5_content:

  data.table. Output of
  [`rhdf5::h5ls()`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html)
  on `f_path`.

- group:

  String. One of `"/obs"` or `"/var"`.

## Value

The number of rows as numeric, `NA` if it cannot be resolved.
