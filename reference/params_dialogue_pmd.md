# Wrapper function for the DIALOGUE decomposition parameters

Stage one of DIALOGUE: the penalised matrix decomposition that turns the
per-cell-type features into multicellular programmes, and the
provisional gene signatures that come off it.

## Usage

``` r
params_dialogue_pmd(
  k = 2L,
  n_permutations = 100L,
  extra_sparse = FALSE,
  abn_c = 15L,
  p_anova = 0.05,
  centre = TRUE,
  cap = 0.01,
  spatial = FALSE,
  n_genes = 200L,
  min_ci = 0.05,
  averaging = c("median", "mean"),
  mcp_assignment_p = 0.1,
  seed = 1234L
)
```

## Arguments

- k:

  Integer. Number of multicellular programmes to extract. Must be at
  least 1.

- n_permutations:

  Integer. Permutations backing the empirical p-value per programme.
  Must be at least 2.

- extra_sparse:

  Boolean. Tune the L1 bound by permutation instead of fixing it at
  `sqrt(p_1) / 2`. Costs ten more fits per permutation.

- abn_c:

  Integer. Minimum cells a sample must contribute, within a cell type,
  before it counts towards the feature-level ANOVA.

- p_anova:

  Numeric. BH-adjusted ANOVA cutoff for keeping a feature. Must be in
  `(0, 1]`.

- centre:

  Boolean. Centre and scale the sample-level feature matrix, then
  winsorise it.

- cap:

  Numeric. Winsorising tail fraction applied to each column. Must be in
  `[0, 0.5)`.

- spatial:

  Boolean. Spatial data: skip the ANOVA feature filter entirely. Niches
  are small, so a feature need not vary across them to be real.

- n_genes:

  Integer. Genes taken per programme per direction when building a
  signature.

- min_ci:

  Numeric. Minimum absolute correlation for a gene to enter a signature.
  Must be in `[0, 1]`.

- averaging:

  String. One of `c("median", "mean")`. How cell-level features are
  collapsed per sample.

- mcp_assignment_p:

  Numeric. Empirical p below which a cell type pair counts as connected
  when deciding which cell types a programme spans. Must be in `(0, 1]`.

- seed:

  Integer. Seed for the permutation null.

## Value

A list with the stage one DIALOGUE parameters.

## Details

The defaults follow upstream's `DLG.get.param`. Two knobs are worth
thinking about before anything else. `k` is how many programmes you are
asking for, and there is no sweep to help you pick it. `n_permutations`
sets the resolution of the empirical p-value: with the default of `100`
the smallest p you can observe is `0.01`, so lower it for a quick look
and leave it alone for anything you intend to believe.

`averaging` is exposed and honoured here. Upstream takes the same
argument and then ignores it, hard-coding column medians, so `"median"`
is what every published DIALOGUE run actually used.

## References

Jerby-Arnon & Regev, Nature Biotechnology, 2022
