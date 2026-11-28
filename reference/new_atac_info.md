# Generate a new `AtacInfo` object

What an ATAC store holds and how far it was processed. The kind, count
mode and TF-IDF flag are read from the store; the LSI fields come from
[`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md)
and are saved with the object.

## Usage

``` r
new_atac_info(kind, count_mode, tf_idf = FALSE)
```

## Arguments

- kind:

  String. One of `c("peak", "tile")`.

- count_mode:

  String. One of `c("insertion", "paired-insertion", "fragment")`. What
  one count is: up to two insertions per fragment, one per fragment and
  feature for the latter two.

- tf_idf:

  Boolean. Whether the store is TF-IDF normalised.

## Value

An `AtacInfo` object, a list with

- kind - `"peak"` or `"tile"`.

- count_mode - What one count is.

- tf_idf - Whether the store is TF-IDF normalised.

- depth_cor - Correlation of each LSI component with the library size,
  all components including any dropped one. `NULL` before the LSI.

- lsi_dropped - 0-based LSI components left out of the `"lsi"`
  embedding.

- lsi_features - 0-based feature rows the LSI was fitted on.

- lsi_loadings - Features x components loadings of the LSI.

- lsi_singular_values - Singular values of the LSI.
