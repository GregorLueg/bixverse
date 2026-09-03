# Get the fitted model

**\[deprecated\]** The neighbourhood test moved into Rust, so there is
no `DGEGLM` to hand back any more. Everything the fit was consulted for
now sits in the results table, see
[`get_differential_abundance_res()`](https://gregorlueg.github.io/bixverse/reference/get_differential_abundance_res.md).

## Usage

``` r
get_model_fit(x)

# S3 method for class 'miloR'
get_model_fit(x)
```

## Arguments

- x:

  An object from which to get the fitted model from.

## Value

`NULL`, with a deprecation warning.
