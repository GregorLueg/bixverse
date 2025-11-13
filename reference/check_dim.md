# Helper to assert same number of samples

Helper to assert same number of samples

## Usage

``` r
check_dim(object, no_samples)
```

## Arguments

- object:

  `snf` class.

- no_samples:

  Integer. New number of samples

## Value

The object. If no dimensions were found in the object yet, it will be
added. If the addition of the new data has dimensions that do not fit
with what has been stored in the class so far, it will throw an error.
