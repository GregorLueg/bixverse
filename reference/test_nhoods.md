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
#>    Nhood       logFC   logCPM            F    PValue       FDR SpatialFDR
#>    <int>       <num>    <num>        <num>     <num>     <num>      <num>
#> 1:     1  0.42390689 13.85440 0.3958957378 0.5296492 0.9455350  0.9459114
#> 2:     2  0.11686896 14.17130 0.0458927514 0.8305024 0.9455350  0.9459114
#> 3:     3 -0.42239717 13.76630 0.3524301245 0.5531450 0.9455350  0.9459114
#> 4:     4 -0.37967351 13.85711 0.3177198405 0.5733617 0.9455350  0.9459114
#> 5:     5  0.72784199 13.61784 1.0455122319 0.3078362 0.9455350  0.9459114
#> 6:     6  0.02071073 13.71601 0.0007936764 0.9775417 0.9775417  0.9775417

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
