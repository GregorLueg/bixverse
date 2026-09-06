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

  Named nested list. Every element must itself be a list with at least a
  `"pos"` element holding the gene identifiers, and optionally a `"neg"`
  one. A bare character vector is not accepted. The gene identifiers
  need to be part of the variables of the object.

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

## Examples

``` r
# a signed signature alongside a plain one
sc <- demo_single_cells()
gs_list <- list(
  programme_a = list(
    pos = get_gene_names(sc)[1:10],
    neg = get_gene_names(sc)[11:20]
  ),
  programme_b = list(pos = get_gene_names(sc)[21:30])
)
res <- vision_sc(sc, gs_list = gs_list, .verbose = FALSE)
dim(res)
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
