# Run NEBULA on single cells

Fits NEBULA's negative binomial gamma mixed model, which is the test to
reach for when the cells are not independent because they came from a
handful of donors. The variance is split into a subject-level random
effect and a cell-level overdispersion, so the donor structure is
modelled rather than ignored, and the nominal false discovery rate
holds.

The design is a formula evaluated against the obs table, so anything in
there can go into the model. Cells with a missing design or subject
value are dropped with a warning.

Genes are streamed and fitted in batches, which is exact: NEBULA is
gene-independent. A fit is milliseconds to seconds per gene though, so
running it over the full gene axis is rarely what you want. Restrict
`genes_to_use` to the highly variable genes or a candidate set.

## Usage

``` r
nebula_sc(
  object,
  subject_col,
  design,
  coef = NULL,
  contrast = NULL,
  genes_to_use = NULL,
  offset = NULL,
  nebula_params = params_nebula(),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- subject_col:

  String. The column in the obs table holding the subject (donor)
  identifier. This is what the random effect is over. Not the same thing
  as a sample or a batch, unless they happen to coincide.

- design:

  Formula. The experimental design, evaluated against the obs table,
  e.g. `~ condition` or `~ condition + age`. Include the intercept.

- coef:

  Optional integer or character. Which coefficient of the design the
  Wald test reports, as a 1-based column position or a column name.
  Defaults to the last column.

- contrast:

  Optional numeric vector. One weight per design column. Mutually
  exclusive with `coef`.

- genes_to_use:

  Optional character vector. The genes to fit. Defaults to every gene in
  the object, which is usually too many.

- offset:

  Optional numeric vector. Strictly positive scaling factor per cell,
  aligned to the cells that survive the design. Defaults to `NULL`,
  which uses the library sizes.

- nebula_params:

  A list, see
  [`params_nebula()`](https://gregorlueg.github.io/bixverse/reference/params_nebula.md).
  The list has the following parameters:

  - nebula_method - String. One of `c("ln", "hl")`.

  - min_sigma, max_sigma - Numeric. Bounds on the subject-level
    overdispersion.

  - min_phi, max_phi - Numeric. Bounds on the cell-level overdispersion.

  - cutoff_cell - Numeric. When to refit both overdispersions.

  - kappa - Numeric. When to trust the stage-one subject overdispersion.

  - cpc - Numeric. Minimum mean count per cell for a gene to be tested.

  - mincp - Integer. Minimum number of cells expressing a gene.

  - reml - Boolean. Restricted maximum likelihood.

  - eps - Numeric. Optimiser stopping tolerance.

  - gene_batch_size - Integer. Genes read and fitted per batch.

  - shrink_dispersion - Boolean. Empirical Bayes shrinkage of the
    cell-level overdispersions.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `ScNebula` class, see
[`new_sc_nebula_res()`](https://gregorlueg.github.io/bixverse/reference/new_sc_nebula_res.md),
with

- results - data.table. One row per gene that survived NEBULA's
  expression filter, with the Wald test and both overdispersions.

- coefficients - Numeric matrix of genes x coefficients.

- se - Numeric matrix of genes x coefficients.

- params - List. The parameters the run used.

## References

He, et al., Commun Biol, 2021
