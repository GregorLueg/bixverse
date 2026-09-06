# Calculate module activity scores

Implements the approach from Tirosh et al into bixverse, i.e., the
AddModuleScore functionality from Seurat. For each module (gene set),
computes the average expression of genes in the set minus the average
expression of randomly selected control genes from the same expression
bins. Genes are binned based on their average expression across cells to
ensure controls are expression-matched.

## Usage

``` r
module_scores_sc(
  object,
  gs_list,
  n_bins = 24L,
  n_ctrl = 100L,
  seed = 42L,
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- gs_list:

  Named list. The elements have the gene identifiers of the respective
  gene sets.

- n_bins:

  Integer. Number of bins to use. Defaults to `24L`.

- n_ctrl:

  Integer. Number of control genes to use per gene for each gene set.
  Defaults to `100L`.

- seed:

  Integer. The random seed.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

Returns a `ScMatrixRes` with the module scores.

## References

Tirosh et al, Science (2016)

## Examples

``` r
# score the three planted marker programmes per cell
sc <- demo_single_cells()
gs_list <- list(
  type_1 = sprintf("gene_%02d", 1:10),
  type_2 = sprintf("gene_%02d", 11:20),
  type_3 = sprintf("gene_%02d", 21:30)
)
res <- module_scores_sc(sc, gs_list = gs_list, n_ctrl = 5L, .verbose = FALSE)
dim(res)
#> [1] 500   3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
