# Run NEBULA on meta cells

Fits NEBULA's negative binomial gamma mixed model over aggregated
counts. The arithmetic is
[`nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/nebula_sc.md)
verbatim, only the counts come out of memory rather than the streamed
store, and the fit is much cheaper because there are far fewer rows.

What changes is the interpretation, not the numerics. The cell-level
overdispersion becomes the spread between aggregates within a subject
rather than between cells, so it is smaller and it absorbs whatever the
aggregation smoothed away. The subject-level term keeps its meaning.
Read the two as a variance decomposition over meta cells and do not
compare them against a single cell run.

## Usage

``` r
nebula_mc(
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

  `MetaCells` class.

- subject_col:

  String. The column in the obs table holding the subject (donor)
  identifier. This is what the random effect is over.

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
  the object.

- offset:

  Optional numeric vector. Strictly positive scaling factor per meta
  cell, aligned to the meta cells that survive the design. Defaults to
  `NULL`, which uses the aggregated library sizes.

- nebula_params:

  A list, see
  [`params_nebula()`](https://gregorlueg.github.io/bixverse/reference/params_nebula.md).
  See
  [`nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/nebula_sc.md)
  for the individual elements.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `ScNebula` class, see
[`new_sc_nebula_res()`](https://gregorlueg.github.io/bixverse/reference/new_sc_nebula_res.md).

## References

He, et al., Commun Biol, 2021
