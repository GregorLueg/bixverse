# Apply a pipeline to a single cell object

Runs each step in order. The first step receives `object`; subsequent
steps receive the result of the previous step. Errors propagate; nothing
is caught.

## Usage

``` r
apply_pipeline(pipeline, object)
```

## Arguments

- pipeline:

  `ScPipeline`.

- object:

  `SingleCells` or `SingleCellsSubset`. Dispatch happens inside each
  step's underlying generic, so the same pipeline works on either.

## Value

The object after all steps have run.
