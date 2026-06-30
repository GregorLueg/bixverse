# Pathway activity functions

## Single-sample pathway activity

A common question in transcriptomics is whether the genes belonging to a
given pathway are concordantly up- or down-regulated in each sample.
Several algorithms tackle this by collapsing a gene-by-sample expression
matrix into a pathway-by-sample score matrix, the most widely used being
[GSVA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)
and [ssGSEA](https://www.nature.com/articles/nature08460). `biverse`
provides fast versions of both with some tricks under the hood to make
them faster than the original implementations. Moreover, if you have
multiple contrasts (different drugs vs. DMSO for example). `bixverse`
also provides Rust-accelerated implementations from
[mitch](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06856-9).

``` r

library(bixverse)
```

## Data

Following the GSVA vignette convention, we generate some synthetic
expression data (Gaussian-distributed to mimic normalised microarray or
RNA-seq log counts, and Poisson-distributed to mimic raw counts)
together with a collection of random gene sets.

``` r

p <- 10000L # genes
n <- 100L # samples

# Gaussian expression matrix
X <- matrix(
  rnorm(p * n),
  nrow = p,
  dimnames = list(paste0("g", seq_len(p)), paste0("s", seq_len(n)))
)

# Poisson count matrix
X_counts <- matrix(
  rpois(p * n, lambda = 10),
  nrow = p,
  dimnames = list(paste0("g", seq_len(p)), paste0("s", seq_len(n)))
)
storage.mode(X_counts) <- "numeric"

# Random gene sets of varying size
gs <- as.list(sample(10:100, size = 250, replace = TRUE))
gs <- lapply(
  gs,
  function(n, p) paste0("g", sample(seq_len(p), size = n, replace = FALSE)),
  p
)
names(gs) <- paste0("gs", seq_along(gs))
```

## GSVA

### Gaussian kernel

The Gaussian kernel version of GSVA is appropriate for continuous
expression values (log-CPM, microarray intensities, etc.). Running it in
`bixverse` is a single call:

``` r

bixverse_res_gaussian <- calc_gsva(
  exp = X,
  pathways = gs,
  kernel = "gaussian"
)

bixverse_res_gaussian[1:5, 1:5]
#>              s1          s2          s3          s4          s5
#> gs1 -0.02762868 -0.06164412 -0.22442119 -0.18778701  0.02012806
#> gs2  0.26850160  0.02513455  0.04308302  0.08502409 -0.23815456
#> gs3  0.39306189 -0.08191868 -0.22790010  0.47293974  0.11599080
#> gs4  0.11971514 -0.22561081 -0.14926845 -0.10618037  0.08759237
#> gs5 -0.11585761  0.24794621  0.15359424  0.16720155  0.24310646
```

If the original `GSVA` Bioconductor package is available we can verify
that the results are essentially identical. Minor differences in
numerical precision arise from optimisations on the Rust side, but
per-pathway correlations between the two implementations are
consistently above 0.99.

``` r

library(GSVA)

gsvaPar <- gsvaParam(X, gs)
gsva_res_gaussian <- as.matrix(gsva(gsvaPar, verbose = FALSE))

correlations <- diag(cor(bixverse_res_gaussian, gsva_res_gaussian))

print(sprintf(
  "All pathway correlations >= 0.99: %s (min: %.4f)",
  all(correlations >= 0.99),
  min(correlations)
))
#> [1] "All pathway correlations >= 0.99: TRUE (min: 1.0000)"
```

Let’s compare speed (differences will be more pronounced, the more
computational fire power your system has.)

``` r

microbenchmark::microbenchmark(
  gsva = {
    gsvaPar <- gsvaParam(X, gs)
    gsva(gsvaPar, verbose = FALSE)
  },
  bixverse = calc_gsva(exp = X, pathways = gs, kernel = "gaussian"),
  times = 3L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 2511.8373 2519.1374 2530.0666 2526.4375 2539.1813 2551.9250     3
#>  bixverse  555.9847  556.5191  557.0972  557.0535  557.6534  558.2534     3
```

### Poisson kernel

For count data a Poisson kernel is more appropriate. The only change is
setting `kernel = "poisson` (to note: prior to `"0.4.1"`, this was
controlled via a Boolean flag `gaussian = FALSE` - this still works, but
will throw a deprecation warning now):

``` r

bixverse_res_poisson <- calc_gsva(
  exp = X_counts,
  pathways = gs,
  kernel = "poisson"
)

bixverse_res_poisson[1:5, 1:5]
#>              s1           s2          s3          s4           s5
#> gs1  0.17386353 -0.127472014  0.07643971 -0.18925777  0.009891965
#> gs2 -0.04066373 -0.091788267  0.15776229 -0.14526685  0.068958105
#> gs3 -0.01870068 -0.003060888 -0.27571781  0.08929506 -0.100720476
#> gs4 -0.14410134 -0.203525162  0.04604448  0.13324940 -0.023074976
#> gs5  0.04771749  0.216379673 -0.07698620 -0.11687178  0.040224566
```

And vs. the original

``` r

gsvaParPoisson <- gsvaParam(X_counts, gs, kcdf = "Poisson")
gsva_res_poisson <- as.matrix(gsva(gsvaParPoisson, verbose = FALSE))

correlations <- diag(cor(bixverse_res_poisson, gsva_res_poisson))

print(sprintf(
  "All pathway correlations >= 0.99: %s (min: %.4f)",
  all(correlations >= 0.99),
  min(correlations)
))
#> [1] "All pathway correlations >= 0.99: TRUE (min: 1.0000)"
```

And speed:

``` r

microbenchmark::microbenchmark(
  gsva = {
    gsvaPar <- gsvaParam(X_counts, gs, kcdf = "Poisson")
    gsva(gsvaPar, verbose = FALSE)
  },
  bixverse = calc_gsva(exp = X_counts, pathways = gs, kernel = "poisson"),
  times = 3L
)
#> Unit: seconds
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 18.898024 18.898101 18.912399 18.898178 18.919586 18.940994     3
#>  bixverse  2.684689  2.685066  2.685917  2.685443  2.686532  2.687621     3
```

## ssGSEA

ssGSEA, first described in [Barbie et
al.](https://www.nature.com/articles/nature08460), takes a similar
approach but does not apply a kernel-based normalisation across samples.
`bixverse` provides a Rust-optimised version here as well:

``` r

bixverse_res_ssgsea <- calc_ssgsea(
  exp = X,
  pathways = gs
)

bixverse_res_ssgsea[1:5, 1:5]
#>            s1          s2          s3          s4          s5
#> gs1 0.1206507  0.01275135  0.01902510 -0.03130446  0.13723916
#> gs2 0.3031894  0.12862560  0.08015432  0.17991118 -0.04373609
#> gs3 0.2498230  0.20226982 -0.03227660  0.38983382  0.18276653
#> gs4 0.1466774 -0.01969500  0.02983841  0.03402698  0.16856867
#> gs5 0.0503750  0.21660765  0.25548401  0.12477599  0.24650876
```

Let’s compare again against the GSVA version:

``` r

ssgseaPar <- ssgseaParam(X, gs)
ssgsea_res <- as.matrix(gsva(ssgseaPar, verbose = FALSE))

correlations <- diag(cor(bixverse_res_ssgsea, ssgsea_res))

print(sprintf(
  "All pathway correlations >= 0.99: %s (min: %.4f)",
  all(correlations >= 0.99),
  min(correlations)
))
#> [1] "All pathway correlations >= 0.99: TRUE (min: 1.0000)"
```

And check the underlying speed differences:

``` r

microbenchmark::microbenchmark(
  gsva = {
    ssgseaPar <- ssgseaParam(X, gs)
    gsva(ssgseaPar, verbose = FALSE)
  },
  bixverse = calc_ssgsea(exp = X, pathways = gs),
  times = 3L
)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>      gsva 606.8413 608.1971 612.3263 609.5530 615.0688 620.5846     3
#>  bixverse 121.6192 121.6737 121.7655 121.7282 121.8387 121.9491     3
```

## singscore

[singscore](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2435-4)
takes a different angle than GSVA or ssGSEA: rather than estimating a
null distribution over samples or genes, it scores each sample
independently using only the ranks of the genes in the signature. This
makes it well-suited to scoring single samples in isolation and to
working with directional signatures (separate up- and down-regulated
gene sets). `bixverse` provides a Rust-accelerated implementation.

### Ranking

singscore operates on a ranked expression matrix rather than the raw
values, so the first step is to compute per-sample ranks:

``` r

ranks <- calc_singscore_rank(exp = X)

ranks[1:5, 1:5]
#>      s1   s2   s3   s4   s5
#> g1 1855 2835 9128 4597 2852
#> g2 4384 7165 4976 7519 7225
#> g3 7194 6390 9596 4885  402
#> g4 5013 1200 6193 8261  508
#> g5 2584 3629 7813 7062  950
```

If your data comes from a small targeted panel (NanoString, RT-qPCR)
rather than a transcriptome-wide assay, you can pass a set of stable
genes to `calc_singscore_rank` to use the stable-gene ranking approach
from [Bhuva et
al.](https://academic.oup.com/nar/article/48/16/e97/5882020):

``` r

ranks_stable <- calc_singscore_rank(
  exp = X,
  stable_genes = c("g1", "g2", "g3", "g4", "g5")
)
```

The `stable` attribute on the rank matrix is read automatically by the
downstream scoring functions to pick the appropriate bounds formula.

### Scoring a single signature

For a single signature, `calc_singscore` returns one score and
dispersion per sample. Below we use one of our random gene sets as the
up-regulated set and another as the down-regulated set:

``` r

up_set <- gs[[1]]
down_set <- gs[[2]]

singscore_res <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set
)

head(singscore_res)
#>    total_score total_dispersion    up_score up_dispersion  down_score
#>          <num>            <num>       <num>         <num>       <num>
#> 1: -0.09600291         3191.672  0.00111725      3647.201 -0.09712016
#> 2: -0.05100830         3439.266 -0.03969083      2957.791 -0.01131747
#> 3: -0.01673365         3851.801 -0.03848764      3718.366  0.02175399
#> 4: -0.10434697         3813.994 -0.07372398      3328.442 -0.03062299
#> 5:  0.08950395         3791.385  0.01082092      3900.726  0.07868304
#> 6:  0.05113496         2256.521  0.01416039      2982.996  0.03697457
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        2736.142        s1
#> 2:        3920.742        s2
#> 3:        3985.235        s3
#> 4:        4299.546        s4
#> 5:        3682.043        s5
#> 6:        1530.045        s6
```

If the direction of the signature is unknown (a gene ontology term, for
example), pass only `up_set` and set `known_direction = FALSE`. The
scoring function then transforms ranks around their median so that genes
at either extreme contribute symmetrically.

### Permutation testing

To assess whether an observed score is larger than would be expected by
chance for a random gene set of the same size, set `n_permutations` to a
positive integer. The function draws random gene sets matching the size
of the real signature, scores them, and returns empirical one-tailed
p-values alongside the scores. The full null distribution is attached as
an attribute for inspection or plotting.

``` r

singscore_perm <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set,
  n_permutations = 1000L,
  seed = 42L
)

head(singscore_perm)
#>    total_score total_dispersion    up_score up_dispersion  down_score
#>          <num>            <num>       <num>         <num>       <num>
#> 1: -0.09600291         3191.672  0.00111725      3647.201 -0.09712016
#> 2: -0.05100830         3439.266 -0.03969083      2957.791 -0.01131747
#> 3: -0.01673365         3851.801 -0.03848764      3718.366  0.02175399
#> 4: -0.10434697         3813.994 -0.07372398      3328.442 -0.03062299
#> 5:  0.08950395         3791.385  0.01082092      3900.726  0.07868304
#> 6:  0.05113496         2256.521  0.01416039      2982.996  0.03697457
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        2736.142        s1 0.001
#> 2:        3920.742        s2 0.001
#> 3:        3985.235        s3 1.000
#> 4:        4299.546        s4 0.004
#> 5:        3682.043        s5 0.001
#> 6:        1530.045        s6 1.000

dim(attr(singscore_perm, "null_distribution"))
#> [1] 1000  100
```

The seed argument makes the permutations reproducible.

### Scoring many signatures

When you have many signatures to score at once, `calc_singscore_multi`
avoids the per-call overhead and parallelises across gene sets. It
accepts a named list of up-regulated sets, optionally paired with a list
of down-regulated sets keyed by the same names:

``` r

singscore_multi_res <- calc_singscore_multi(
  ranks = ranks,
  up_pathways = gs[1:10]
)

singscore_multi_res$scores[1:5, 1:5]
#>              s1          s2          s3          s4          s5
#> gs1  0.00111725 -0.03969083 -0.03848764 -0.07372398  0.01082092
#> gs2  0.09712016  0.01131747 -0.02175399  0.03062299 -0.07868304
#> gs3  0.09098645  0.02429491 -0.08111195  0.16290192  0.05431884
#> gs4  0.02223504 -0.06647860 -0.04118359 -0.03504746  0.02920010
#> gs5 -0.02860509  0.06442324  0.06095639  0.02942167  0.08596570
```

The returned list contains a `scores` matrix and a matching
`dispersions` matrix, both with shape gene sets × samples.

### Comparison with singscore

When the original `singscore` package is available we can confirm that
the two implementations agree closely. One detail worth noting:
`singscore` uses `ties.method = "min"` when ranking, while `bixverse`
uses average ranks. For continuous data the difference is negligible,
but it can introduce small discrepancies on heavily tied data.

``` r

library(singscore)

sing_ranks <- rankGenes(X)
sing_res <- simpleScore(
  rankData = sing_ranks,
  upSet = up_set,
  downSet = down_set
)

cor(singscore_res$total_score, sing_res$TotalScore)
```

And the speed difference:

``` r

microbenchmark::microbenchmark(
  singscore = {
    sing_ranks <- rankGenes(X)
    simpleScore(rankData = sing_ranks, upSet = up_set, downSet = down_set)
  },
  bixverse = {
    ranks <- calc_singscore_rank(exp = X)
    calc_singscore(ranks = ranks, up_set = up_set, down_set = down_set)
  },
  times = 5L
)
```

## Multi-contrast enrichment (mitch)

When you have multiple contrasts - say, differential expression results
from several comparisons or different omics layers - you often want to
know which pathways show coordinated changes *across* contrasts. The
[mitch](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06856-9)
method does exactly this: it ranks genes within each contrast, computes
mean ranks per pathway, and then uses a MANOVA test to identify pathways
that are enriched in one or more contrasts simultaneously. `bixverse`
re-implements the core algorithm so that it slots into the same workflow
as the other pathway activity functions… Also, heavily multi-threaded
and designed to go **brrrrrr** in terms of speed.

`calc_mitch` expects a numeric matrix of contrast statistics (genes in
rows, contrasts in columns) and a named list of gene sets. It returns a
`data.table` containing per-pathway MANOVA p-values, FDR-adjusted
p-values, and the individual contrast-level enrichment scores and
p-values.

``` r

set.seed(42L)

contrast_data <- matrix(rnorm(3 * 26), nrow = 26)
colnames(contrast_data) <- sprintf("contrast_%i", 1:3)
rownames(contrast_data) <- letters

gene_sets <- list(
  pathway_A = sample(letters, 4),
  pathway_B = sample(letters, 5),
  pathway_C = sample(letters, 6),
  pathway_D = sample(letters, 7)
)

res <- calc_mitch(
  contrast_mat = contrast_data,
  gene_set_list = gene_sets
)

res
```

Note that pathways smaller than the internal minimum set size are
dropped from the output, and the results are sorted by significance. The
function also validates its input: passing a matrix containing `NA`
values will throw an error rather than silently producing nonsense.

### Comparison with the mitch package

When the original `mitch` package is available we can verify that the
two implementations produce identical results. Here we use the example
data bundled with `mitch`:

``` r

data(myImportedData, genesetsExample, package = "mitch")

mitch_res <- suppressMessages(mitch::mitch_calc(
  myImportedData,
  genesetsExample,
  priority = "significance",
  minsetsize = 5,
  cores = 2
))

bixverse_res <- calc_mitch(
  contrast_mat = as.matrix(myImportedData),
  gene_set_list = genesetsExample
)

# Pathway ordering and FDR values should match exactly
all(bixverse_res$pathway_names == mitch_res$enrichment_result$set) &&
  all.equal(bixverse_res$manova_fdr, mitch_res$enrichment_result$p.adjustMANOVA)
#> [1] TRUE
```

And let’s check the speed differences:

``` r

microbenchmark::microbenchmark(
  mitch = suppressMessages(mitch::mitch_calc(
    myImportedData,
    genesetsExample,
    priority = "significance",
    minsetsize = 5,
    cores = 2
  )),
  bixverse = calc_mitch(
    contrast_mat = as.matrix(myImportedData),
    gene_set_list = genesetsExample
  ),
  times = 5L
)
#> Unit: milliseconds
#>      expr        min         lq       mean     median         uq       max
#>     mitch 143.909298 144.099963 146.126884 145.174249 148.382494 149.06842
#>  bixverse   3.428616   3.683632   4.409473   3.929542   5.163463   5.84211
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
