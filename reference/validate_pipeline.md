# Check that a pipeline can run on a given class

Walks the pipeline and tracks which class each step would receive. Steps
declare the classes their generic has methods for, and most steps hand
back what they were given;
[`step_metacells_sc()`](https://gregorlueg.github.io/bixverse/reference/step_metacells_sc.md)
does not, it turns a `SingleCells`/`SingleCellsSubset` into a
`MetaCells`. Without this you would only find out about a mismatch when
S7 fails to dispatch, i.e. after the expensive steps already ran.

## Usage

``` r
validate_pipeline(pipeline, class)
```

## Arguments

- pipeline:

  `ScPipeline`.

- class:

  String. Class of the object the pipeline would start on. One of
  `c("SingleCells", "SingleCellsSubset", "MetaCells")`.

## Value

Invisibly, the class the pipeline would return.
