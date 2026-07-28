# Append a step to a pipeline

`%>>%` chains pipeline steps. Either side can be a `ScStep` or a
`ScPipeline`; the result is always a `ScPipeline`.

## Usage

``` r
lhs %>>% rhs
```

## Arguments

- lhs:

  `ScPipeline` or `ScStep`.

- rhs:

  `ScStep`.

## Value

A `ScPipeline`.
