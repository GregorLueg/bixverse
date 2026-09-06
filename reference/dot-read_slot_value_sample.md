# Read a sample of values from a matrix slot in an h5ad file

Handles both sparse (group with a `data` dataset) and dense (direct
dataset) storage. Returns NULL if the slot does not exist.

## Usage

``` r
.read_slot_value_sample(f_path, slot_path, h5_content, n_sample = 10000L)
```

## Arguments

- f_path:

  String. Path to the h5ad file.

- slot_path:

  String. Full path of the slot inside the file.

- h5_content:

  data.table. Output of
  [`rhdf5::h5ls()`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html)
  on `f_path`.

- n_sample:

  Integer. Maximum number of values to read.

## Value

Numeric vector of sampled values, or `NULL` if the slot is absent.
