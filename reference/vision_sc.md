# Calculate VISION scores

Calculates an VISION-type scores for pathways based on DeTomaso, et al.
Compared to other score types, you can also calculate delta-type scores
between positive and negative gene indices, think epithelial vs
mesenchymal gene signature, etc.

## Usage

``` r
vision_sc(object, gs_list, streaming = FALSE, .verbose = TRUE)
```

## Arguments

- object:

  `single_cell_exp` class.

- gs_list:

  Named nested list. The elements have the gene identifiers of the
  respective gene sets and have the option to have a `"pos"` and `"neg"`
  gene sets. The names need to be part of the variables of the
  `single_cell_exp` class.

- streaming:

  Boolean. Shall the cell data be streamed in. Useful for larger data
  sets.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Matrix of cells x signatures with the VISION pathway scores as values.

## References

DeTomaso, et al., Nat. Commun., 2019
