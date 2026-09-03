# Fit the NEBULA negative binomial gamma mixed model over meta cells

**\[experimental\]** The arithmetic is the single cell one verbatim,
only the counts come from memory rather than the streamed store. What
changes is the interpretation: the cell-level overdispersion becomes the
spread between aggregates within a subject, not between cells, so it is
smaller and absorbs whatever the aggregation smoothed away. The
subject-level term keeps its meaning. Do not compare the two against a
single cell run.

## Usage

``` r
rs_nebula_mc(
  sparse_data,
  metacells_to_keep,
  gene_indices,
  subject_ids,
  design,
  offset,
  nebula_params,
  verbose
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape (metacells, genes), holding the raw
  counts.

- metacells_to_keep:

  Integer vector. 0-indexed(!) positions of the meta cells to analyse,
  in any order. Must not hold duplicates.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to fit.

- subject_ids:

  Integer vector. 0-indexed(!) subject label per meta cell. One entry
  per row of `sparse_data`, not per element of `metacells_to_keep`.

- design:

  Numeric matrix. Predictors of meta cells x coefficients, rows aligned
  to `metacells_to_keep` and including an intercept.

- offset:

  Optional numeric vector. Strictly positive scaling factor per selected
  meta cell. `NULL` uses the aggregated library sizes.

- nebula_params:

  Named list. The NEBULA parameters, see
  [`params_nebula()`](https://gregorlueg.github.io/bixverse/reference/params_nebula.md),
  plus either `coef` (a 0-indexed(!) coefficient) or `contrast` (one
  weight per coefficient).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the same elements
[`rs_nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nebula_sc.md)
returns.

## References

He, et al., Commun Biol, 2021
