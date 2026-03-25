# Calculate VISION scores (with auto-correlation scores)

Calculates VISION-type scores for pathways based on DeTomaso, et al.
Compared to other score types, you can also calculate delta-type scores
between positive and negative gene indices, think epithelial vs
mesenchymal gene signature, etc. Additionally, this function also
calculates the auto- correlation values, answering the question if a
given signature shows non- random enrichment on the kNN graph. The kNN
graph (and distance measures) will be generated on-the-fly based on the
embedding you wish to use.

## Usage

``` r
vision_w_autocor_sc(
  object,
  gs_list,
  embd_to_use,
  no_embd_to_use = NULL,
  vision_params = params_sc_vision(),
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- gs_list:

  Named nested list. The elements have the gene identifiers of the
  respective gene sets and have the option to have a `"pos"` and `"neg"`
  gene sets. The names need to be part of the variables of the
  `SingleCells` class.

- embd_to_use:

  String. The embedding to use. Whichever you chose, it needs to be part
  of the object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- vision_params:

  List with vision parameters, see
  [`params_sc_vision()`](https://gregorlueg.github.io/bixverse/reference/params_sc_vision.md)
  with the following elements:

  - n_perm - Integer. Number of random permutations

  - n_cluster - Integer. Number of random clusters to generate to
    associate each set with.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- streaming:

  Boolean. Shall the cell data be streamed in. Useful for larger data
  sets.

- random_seed:

  Integer. The random seed.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Matrix of cells x signatures with the VISION pathway scores as values.

## References

DeTomaso, et al., Nat. Commun., 2019
