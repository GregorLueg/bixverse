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
#> 1: gene_01 -0.4345582 0.2976426 -1.4599999 0.1442901 0.9274431
#> 2: gene_02 -0.1984760 0.2907377 -0.6826634 0.4948196 0.9274431
#> 3: gene_03 -0.1901702 0.2949573 -0.6447379 0.5190970 0.9274431
#> 4: gene_04 -0.2830729 0.3538626 -0.7999514 0.4237390 0.9274431
#> 5: gene_05 -0.3958592 0.4129879 -0.9585249 0.3377981 0.9274431
#> 6: gene_06 -0.2923441 0.2664164 -1.0973202 0.2725015 0.9274431
#>    subject_overdispersion cell_overdispersion convergence sigma_at_bound
#>                     <num>               <num>       <int>         <lgcl>
#> 1:             0.05448843            1.150486           1          FALSE
#> 2:             0.05224946            1.077529           1          FALSE
#> 3:             0.05974837            1.021993         -10          FALSE
#> 4:             0.10133330            1.170186           1          FALSE
#> 5:             0.14619244            1.391691           1          FALSE
#> 6:             0.03321586            1.091282           1          FALSE
#>    cell_overdispersion_shrunk
#>                         <num>
#> 1:                   1.119392
#> 2:                   1.063358
#> 3:                   1.001475
#> 4:                   1.108070
#> 5:                   1.310626
#> 6:                   1.066635

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
