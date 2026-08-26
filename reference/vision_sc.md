# Calculate VISION scores

Calculates an VISION-type scores for pathways based on DeTomaso, et al.
Compared to other score types, you can also calculate delta-type scores
between positive and negative gene indices, think epithelial vs
mesenchymal gene signature, etc.

## Usage

``` r
vision_sc(object, gs_list, streaming = NULL, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- gs_list:

  Named nested list. The elements have the gene identifiers of the
  respective gene sets and have the option to have a `"pos"` and `"neg"`
  gene sets. The names need to be part of the variables of the object.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect. Ignored when applied to
  `MetaCells`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The VISION scores in form of a matrix that is cells x gene sets or as
`ScMatrixRes` pending the input.

## References

DeTomaso, et al., Nat. Commun., 2019
