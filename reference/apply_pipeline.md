# Apply a pipeline to a single cell object

Runs each step in order. The first step receives `object`; subsequent
steps receive the result of the previous step. The chain is validated
against the class of `object` before anything runs, see
[`validate_pipeline()`](https://gregorlueg.github.io/bixverse/reference/validate_pipeline.md).
Errors propagate; nothing is caught.

## Usage

``` r
apply_pipeline(pipeline, object)
```

## Arguments

- pipeline:

  `ScPipeline`.

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells`. Dispatch happens
  inside each step's underlying generic, so the same pipeline works on
  any class its steps have methods for.

## Value

The object after all steps have run.
