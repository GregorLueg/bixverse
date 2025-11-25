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
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

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

  Boolean. Shall the cell and gene data be streamed in. Useful for
  larger data sets.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

description

## References

Tirosh et al, Science (2016)
