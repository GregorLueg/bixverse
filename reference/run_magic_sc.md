# Impute a subset of genes with MAGIC

MAGIC smooths counts over the kNN graph: a row-stochastic diffusion
operator is built from the graph stored on the object and applied
`n_steps` times, so every cell becomes a weighted average of its
neighbourhood. It makes sparse marker genes readable on an embedding and
gives gene trends something less speckled to fit. For details, please
refer to van Dijk, et al.

Deliberately restricted to a gene subset. The output is dense and the
binary store is never touched, so pick the genes you want to plot or fit
and impute those. Genome-wide imputation is not offered.

Each call replaces the whole layer. Read it back with
[`get_magic()`](https://gregorlueg.github.io/bixverse/reference/get_magic.md),
or hand it to any of the extractors with `layer = "magic"`.

## Usage

``` r
run_magic_sc(
  object,
  features,
  modality = c("rna", "adt", "wnn"),
  magic_params = params_sc_magic(),
  .verbose = TRUE
)
```

## Arguments

- object:

  One of `SingleCells`, `SingleCellsSubset` or `SingleCellsMultiModal`.
  `MetaCells` are already aggregated and are not supported.

- features:

  Character vector. The genes to impute. Anything not found in the
  object is dropped with a warning.

- modality:

  String. One of `c("rna", "adt", "wnn")`. Which kNN graph does the
  smoothing. Anything but `"rna"` requires a `SingleCellsMultiModal`
  object. The counts themselves always come from the RNA store, and the
  layer is always written to the RNA cache.

- magic_params:

  List. See
  [`params_sc_magic()`](https://gregorlueg.github.io/bixverse/reference/params_sc_magic.md).

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The object with an `ScMagic` layer in its cache.

## What this does to your results

Neighbourhoods overlap, so smoothing over them makes nearby cells
resemble each other in every gene at once. Gene-gene correlation after
imputation is inflated, badly, and the inflation is a property of the
graph rather than of the biology. Do not feed imputed values into
Hotspot, SCENIC, differential correlation or CoReMo: those methods
measure exactly the quantity MAGIC manufactures. Visualisation and trend
fitting are what it is for.

The operator preserves per-cell mass, so imputed values sit on the scale
of whatever went in. Imputing raw counts and imputing log-normalised
counts are different operations rather than the same one rescaled, see
the `layer` element of
[`params_sc_magic()`](https://gregorlueg.github.io/bixverse/reference/params_sc_magic.md).

## References

van Dijk, et al., Cell, 2018.
