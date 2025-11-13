# Helper to generate one-hot encodings

Helper function that generate one-hot encoding with the option of
unlabelled data. (Denoted as `NA` in the labels vector.)

## Usage

``` r
one_hot_encode(labels)
```

## Arguments

- labels:

  String vector. `NA`s signify unlabelled data.

## Value

One-hot encoded matrix.
