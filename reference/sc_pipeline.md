# Construct an empty single cell pipeline

Linear container of `ScStep`s. Append steps with `%>>%` and execute with
[`apply_pipeline()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline.md).
Pipelines are inert until applied; steps can be inspected via
`pipeline$steps`.

## Usage

``` r
sc_pipeline()
```

## Value

An empty `ScPipeline` object.
