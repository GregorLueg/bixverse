# Transform data.tables into matrices for distance calculations

This function can take in data.tables with categorical and continuous
data and will transform them into matrices for subsequent distance
calculations based on Hamming distance (everything categorical) or for
Gower distance (mixed types).

## Usage

``` r
prep_data_gower_hamming_dist(dt, sample_names)
```

## Arguments

- dt:

  data.table. The data.table you want to prepare for Gower distance
  calculation in Rust. Assumes samples x features.

- sample_names:

  String. The sample names. Needs to be same length as `nrow(dt)`.

## Value

List with the following items:

- dat - A numerical matrix ready for Gower distance calculations across
  samples based on mixed features.

- is_cat - Boolean vector, indicating which columns are categorical.
