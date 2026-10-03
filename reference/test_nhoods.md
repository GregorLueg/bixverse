# Test neighbourhoods for differential abundance

Performs differential abundance testing on single-cell neighbourhoods
with edgeR's quasi-likelihood negative binomial framework, implemented
in Rust via the `edge-rs` crate. A generalised linear model is fitted to
the neighbourhood counts, one coefficient or contrast is tested, and the
spatial FDR correction accounts for the fact that neighbourhoods overlap
and their tests are therefore not independent.

`filterByExpr()` is off here and cannot be turned on. It is a gene
expression heuristic and means nothing for a neighbourhood. Use
`min_mean` if you want to drop sparsely populated neighbourhoods.

## Usage

``` r
test_nhoods(
  x,
  design,
  design_df,
  coef = NULL,
  contrast = NULL,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "logMS"),
  min_mean = 0,
  robust = TRUE,
  legacy = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
)

# S3 method for class 'miloR'
test_nhoods(
  x,
  design,
  design_df,
  coef = NULL,
  contrast = NULL,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "logMS"),
  min_mean = 0,
  robust = TRUE,
  legacy = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
)
```

## Arguments

- x:

  `miloR` object for which to run the differential abundance analysis.

- design:

  Formula for the experimental design, e.g. `~ grps`.

- design_df:

  data.frame. The metadata used to build the model matrix. Its rownames
  need to cover the sample names of the neighbourhood counts.

- coef:

  Optional integer or character. Which coefficient(s) of the design to
  drop from the null model, given as 1-based column positions or column
  names. Defaults to the last column, as edgeR does.

- contrast:

  Optional numeric vector or matrix. Weights over the design columns.
  Mutually exclusive with `coef`.

- norm_method:

  String. Library size normalisation. One of
  `c("TMM", "TMMwsp", "RLE", "upperquartile", "logMS")`. Defaults to
  `"TMM"`. `"logMS"` is Milo's own name for leaving every factor at one.

- min_mean:

  Numeric. Minimum mean count across samples. Neighbourhoods below it
  are dropped. Defaults to `0` (no filtering).

- robust:

  Logical. Robust estimation of the quasi-likelihood dispersion.
  Defaults to `TRUE`.

- legacy:

  Logical. Take edgeR's pre-4.0 quasi-likelihood pipeline, which runs
  `estimateDisp()` and applies the Poisson bound. Defaults to `TRUE`, so
  this keeps matching what Milo itself does.

- fdr_weighting:

  String. Spatial FDR weighting scheme. One of
  `c("k-distance", "graph-overlap", "none")`. `"k-distance"` weights by
  the distance to the k-th nearest neighbour, `"graph-overlap"` by the
  number of cells shared with other neighbourhoods. Defaults to
  `"k-distance"`.

## Value

The `miloR` object with the differential abundance results added.

## References

Dann, et al., Nat Biotechnol, 2022; Chen, Lun and Smyth, F1000Research,
2016

## Examples

``` r
# differential abundance of neighbourhoods across two sample groups
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 500L,
    n_genes = 50L,
    n_samples = 6L,
    sample_bias = "even"
  )
)
milo <- get_miloR_abundances_sc(
  sc,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(k_refine = 10L),
  .verbose = FALSE
)
design_df <- data.frame(
  grp = rep(c("a", "b"), each = 3),
  row.names = sprintf("sample_%i", 1:6)
)
milo <- test_nhoods(milo, design = ~grp, design_df = design_df)
head(get_differential_abundance_res(milo))
#>    Nhood       logFC   logCPM           F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>       <num>     <num>     <num>      <num>
#> 1:     1  0.43394889 13.84157 0.433206314 0.5108705 0.9640042  0.9599185
#> 2:     2  0.12629598 14.15793 0.053925002 0.8165104 0.9640042  0.9644670
#> 3:     3 -0.53221669 13.79775 0.599966995 0.4391353 0.9640042  0.9598125
#> 4:     4 -0.36980693 13.84297 0.313858391 0.5756954 0.9640042  0.9599185
#> 5:     5  0.73814666 13.60314 1.138735255 0.2936231 0.9640042  0.9598125
#> 6:     6  0.03387627 13.70417 0.002060917 0.9638176 0.9754298  0.9772304

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
