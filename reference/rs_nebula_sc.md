# Fit the NEBULA negative binomial gamma mixed model over single cells

**\[experimental\]** Fits NEBULA to every requested gene, streaming the
counts out of the gene-major store in batches. NEBULA is
gene-independent, so the batching changes nothing about the answer.
Cells do not have to arrive grouped by subject: the Rust side sorts them
and permutes the design and offsets to match.

## Usage

``` r
rs_nebula_sc(
  f_path_genes,
  f_path_cells,
  cells_to_keep,
  gene_indices,
  subject_ids,
  design,
  offset,
  nebula_params,
  verbose
)
```

## Arguments

- f_path_genes:

  String. Path to the `counts_genes.bin` file.

- f_path_cells:

  String. Path to the `counts_cells.bin` file. Only read when `offset`
  is `NULL`, to take the library sizes.

- cells_to_keep:

  Integer vector. 0-indexed(!) global positions of the cells to analyse,
  in any order. Must not hold duplicates.

- gene_indices:

  Integer vector. 0-indexed(!) positions of the genes to fit.

- subject_ids:

  Integer vector. 0-indexed(!) subject label per global cell. One entry
  per cell in the store, not per cell in `cells_to_keep`.

- design:

  Numeric matrix. Predictors of cells x coefficients, rows aligned to
  `cells_to_keep` and including an intercept.

- offset:

  Optional numeric vector. Strictly positive scaling factor per selected
  cell, aligned to `cells_to_keep`. `NULL` uses the library sizes.

- nebula_params:

  Named list. The NEBULA parameters, see
  [`params_nebula()`](https://gregorlueg.github.io/bixverse/reference/params_nebula.md),
  plus either `coef` (a 0-indexed(!) coefficient) or `contrast` (one
  weight per coefficient).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements

- gene_idx - Integer. 0-indexed positions of the genes that survived
  NEBULA's own expression filter.

- coefficients - Numeric matrix of genes x coefficients. The fixed
  effects on the design scale.

- se - Numeric matrix of genes x coefficients. The standard errors.

- subject_overdispersion - Numeric. NEBULA's `sigma^2`.

- cell_overdispersion - Numeric. NEBULA's `phi^-1`.

- cell_overdispersion_shrunk - Numeric or `NULL`. The cell-level
  overdispersion after empirical Bayes shrinkage, when it was requested.

- convergence - Integer. NEBULA's convergence code. At or below `-20` is
  a likely failure.

- sigma_at_bound - Boolean. Whether the subject-level variance finished
  pinned on its lower bound, i.e. the mixed model collapsed to a plain
  negative binomial.

- log_fc - Numeric. Effect of the tested coefficient or contrast, on the
  natural log scale.

- effect_se - Numeric. Standard error of that effect.

- z - Numeric. The Wald statistic.

- p_values - Numeric. Two-sided p-values.

- fdr - Numeric. Benjamini-Hochberg adjusted p-values.

## References

He, et al., Commun Biol, 2021
