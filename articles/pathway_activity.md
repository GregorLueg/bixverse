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
#>               s1          s2          s3          s4          s5
#> gs1  0.033085838  0.01902984  0.40436854  0.33929020  0.17658850
#> gs2  0.169380835  0.11806933  0.25550687  0.22478973  0.03963821
#> gs3  0.004129519  0.16115670 -0.24521413 -0.04166347 -0.09216345
#> gs4 -0.203716202  0.23757052 -0.27424794  0.20844099 -0.11243070
#> gs5  0.022842770 -0.02727209 -0.09559934  0.23403576 -0.01520081
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
#>      gsva 2194.2867 2199.2728 2212.0485 2204.2590 2220.9294 2237.5998     3
#>  bixverse  475.9807  476.0949  476.2356  476.2091  476.3631  476.5171     3
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
#>              s1          s2          s3          s4          s5
#> gs1 -0.22756639 -0.38285289  0.13230390 -0.03494324 -0.05888921
#> gs2  0.28753661  0.08079601  0.19904746  0.05890301 -0.13835029
#> gs3 -0.02863233  0.03018879 -0.38118921 -0.08271132  0.02210995
#> gs4 -0.01996962 -0.27046455  0.17110917 -0.24308069  0.12532893
#> gs5 -0.06260162 -0.01867168 -0.04783679 -0.21887537  0.02359679
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
#>      gsva 14.963662 14.965130 14.965776 14.966598 14.966833 14.967068     3
#>  bixverse  2.294487  2.295263  2.303137  2.296039  2.307462  2.318885     3
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
#>             s1         s2         s3         s4         s5
#> gs1 0.08015883 0.07596695 0.26790038 0.29591002 0.27823038
#> gs2 0.29083159 0.15229980 0.10235086 0.23642447 0.11011357
#> gs3 0.09412273 0.15203126 0.01090233 0.05694404 0.06659647
#> gs4 0.04478869 0.20548938 0.01798974 0.13198845 0.04515403
#> gs5 0.10508155 0.09151564 0.06447149 0.24224452 0.09629457
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
#>      expr       min       lq     mean    median       uq      max neval
#>      gsva 507.30251 507.9038 517.9219 508.50509 523.2316 537.9581     3
#>  bixverse  96.09757  96.2710 105.1508  96.44444 109.6775 122.9105     3
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
#> g1 7230 7633 7080  740 2041
#> g2 9453 5764 2041 8420 7049
#> g3 9183 6031 5561 3980 2155
#> g4 4528 1947 1609 8378 8611
#> g5 5000 7501 9875  470 9169
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
#> 1: -0.07733740         3548.608  0.00553584      3036.369 -0.08287324
#> 2: -0.03479331         3725.038 -0.01858045      3323.994 -0.01621286
#> 3:  0.08621584         3301.755  0.10717794      3091.226 -0.02096211
#> 4:  0.06449462         4380.348  0.11652130      3396.642 -0.05202668
#> 5:  0.08283606         2242.436  0.10228972      2074.161 -0.01945366
#> 6:  0.02466977         2343.994 -0.02786366      2981.513  0.05253342
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        4060.847        s1
#> 2:        4126.082        s2
#> 3:        3512.285        s3
#> 4:        5364.055        s4
#> 5:        2410.711        s5
#> 6:        1706.475        s6
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
#> 1: -0.07733740         3548.608  0.00553584      3036.369 -0.08287324
#> 2: -0.03479331         3725.038 -0.01858045      3323.994 -0.01621286
#> 3:  0.08621584         3301.755  0.10717794      3091.226 -0.02096211
#> 4:  0.06449462         4380.348  0.11652130      3396.642 -0.05202668
#> 5:  0.08283606         2242.436  0.10228972      2074.161 -0.01945366
#> 6:  0.02466977         2343.994 -0.02786366      2981.513  0.05253342
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        4060.847        s1 0.994
#> 2:        4126.082        s2 0.971
#> 3:        3512.285        s3 0.001
#> 4:        5364.055        s4 0.007
#> 5:        2410.711        s5 1.000
#> 6:        1706.475        s6 0.001

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
#>               s1          s2          s3          s4            s5
#> gs1  0.005535840 -0.01858045  0.10717794  0.11652130  0.1022897243
#> gs2  0.082873237  0.01621286  0.02096211  0.05202668  0.0194536595
#> gs3  0.009503012  0.02854418 -0.06191516 -0.01884538 -0.0189056225
#> gs4 -0.051110164  0.05125635 -0.05241747  0.03277449 -0.0245338813
#> gs5  0.014692770 -0.01305090 -0.02046367  0.08486070  0.0001861018
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
#> [1] 1
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
#> Unit: milliseconds
#>       expr      min       lq     mean   median       uq      max neval
#>  singscore 74.00754 74.12351 74.47907 74.36260 74.60238 75.29929     5
#>   bixverse 15.71730 16.98306 17.80628 17.06735 17.10857 22.15512     5
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
#>      expr        min         lq       mean     median         uq        max
#>     mitch 138.847041 141.438128 143.512382 145.566736 145.631903 146.078103
#>  bixverse   2.685146   2.801128   3.830627   2.998903   4.778461   5.889495
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
