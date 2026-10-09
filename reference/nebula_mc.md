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
  cells_to_use = NULL,
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

- cells_to_use:

  Optional character vector. Meta cell identifiers (`meta_cell_id`) to
  fit. Defaults to every meta cell. Identifiers that cannot be matched
  are dropped with a warning.

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

## Examples

``` r
# differential expression over meta cells with a donor random effect
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 2000L,
    n_genes = 50L,
    n_samples = 6L,
    sample_bias = "even"
  )
)
mc <- generate_bt_meta_cells_sc(
  sc,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 100L),
  .verbose = FALSE
)
# the meta cell obs carries no donor labels, so vote them over
sample_id <- sc[["sample_id"]]$sample_id
mc[["sample_id"]] <- vapply(
  mc[[]]$original_cell_idx,
  function(idx) names(which.max(table(sample_id[idx]))),
  character(1)
)
mc[["condition"]] <- ifelse(
  mc[["sample_id"]]$sample_id %in% c("sample_1", "sample_2", "sample_3"),
  "ctr",
  "trt"
)
res <- nebula_mc(
  mc,
  subject_col = "sample_id",
  design = ~condition,
  .verbose = FALSE
)
head(res$results)
#>    gene_id     log_fc effect_se          z   p_value       fdr
#>     <char>      <num>     <num>      <num>     <num>     <num>
#> 1: gene_01 -0.4345586 0.2976433 -1.4599981 0.1442906 0.9274441
#> 2: gene_02 -0.1984758 0.2907375 -0.6826634 0.4948196 0.9274441
#> 3: gene_03 -0.1901700 0.2949566 -0.6447390 0.5190963 0.9274441
#> 4: gene_04 -0.2830725 0.3538614 -0.7999531 0.4237380 0.9274441
#> 5: gene_05 -0.3958592 0.4129879 -0.9585248 0.3377982 0.9274441
#> 6: gene_06 -0.2923441 0.2664166 -1.0973195 0.2725017 0.9274441
#>    subject_overdispersion cell_overdispersion convergence sigma_at_bound
#>                     <num>               <num>       <int>         <lgcl>
#> 1:             0.05448892            1.150487         -10          FALSE
#> 2:             0.05224929            1.077529           1          FALSE
#> 3:             0.05974774            1.021994           1          FALSE
#> 4:             0.10133211            1.170187           1          FALSE
#> 5:             0.14619244            1.391692           1          FALSE
#> 6:             0.03321597            1.091283           1          FALSE
#>    cell_overdispersion_shrunk
#>                         <num>
#> 1:                  0.9756426
#> 2:                  1.0617131
#> 3:                  1.0192770
#> 4:                  1.1052133
#> 5:                  1.3277607
#> 6:                  1.0653603

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
