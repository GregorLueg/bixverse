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
#>             s1           s2          s3           s4          s5
#> gs1 0.08212298  0.100506350  0.33020082 -0.252741209  0.25715813
#> gs2 0.06587374 -0.182273810 -0.03292728  0.002877405 -0.09753566
#> gs3 0.08642498  0.051235777  0.13749580 -0.083741704  0.11711146
#> gs4 0.15038391  0.272433852  0.06685172  0.127011669  0.16671548
#> gs5 0.12388193  0.002305281  0.06007909 -0.120604485 -0.16981973
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
#>      expr       min        lq      mean   median        uq       max neval
#>      gsva 2836.0882 2838.7096 2846.4867 2841.331 2851.6859 2862.0408     3
#>  bixverse  583.4857  583.7248  584.5996  583.964  585.1565  586.3491     3
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
#>              s1          s2          s3           s4          s5
#> gs1  0.14816902 -0.07755545 -0.06370789 -0.056465417  0.02146155
#> gs2  0.02784901 -0.09784672 -0.04150210 -0.019320910  0.04775480
#> gs3  0.15241034 -0.18270173  0.09876472 -0.021430384  0.09842161
#> gs4  0.10767802  0.41797524 -0.11922181  0.306467131  0.14406874
#> gs5 -0.21094428 -0.01813290 -0.10397481 -0.008051912 -0.07532410
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
#>      gsva 19.293961 19.302678 19.308654 19.311396 19.316001 19.320606     3
#>  bixverse  2.869695  2.870278  2.876189  2.870861  2.879436  2.888012     3
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
#>            s1          s2         s3          s4          s5
#> gs1 0.1245710  0.09954493 0.25775115 -0.03874579 0.179849210
#> gs2 0.1511253 -0.01492187 0.10005123  0.09062564 0.035609564
#> gs3 0.1549151  0.09712388 0.15192389  0.03680863 0.114994004
#> gs4 0.1616648  0.29513648 0.09751194  0.09892486 0.147949974
#> gs5 0.1550947  0.11495952 0.09286364  0.03996919 0.001980941
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
#>      gsva 687.6401 698.9856 706.7611 710.3312 716.3216 722.3120     3
#>  bixverse 128.9757 130.0709 136.0039 131.1662 139.5180 147.8699     3
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
#> g1 3721 9647 9941 8417  995
#> g2 8775 4311 4354 5923  129
#> g3 7454 9063 8658 1636 5623
#> g4  592 1487 6705 9339 1961
#> g5 1668  364 3483 9751  434
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
#>    total_score total_dispersion     up_score up_dispersion   down_score
#>          <num>            <num>        <num>         <num>        <num>
#> 1: -0.01050262         3552.686  0.021100330      3527.111 -0.031602951
#> 2:  0.06882734         2927.027  0.013989979      2894.040  0.054837357
#> 3:  0.10570442         3258.389  0.102853985      2748.745  0.002850436
#> 4: -0.06486347         3031.551 -0.062987209      2492.254 -0.001876258
#> 5:  0.08746122         3580.484  0.049451496      3851.801  0.038009725
#> 6: -0.01276267         3915.182  0.009125326      4338.094 -0.021887995
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        3578.260        s1
#> 2:        2960.015        s2
#> 3:        3768.034        s3
#> 4:        3570.847        s4
#> 5:        3309.168        s5
#> 6:        3492.270        s6
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
#>    total_score total_dispersion     up_score up_dispersion   down_score
#>          <num>            <num>        <num>         <num>        <num>
#> 1: -0.01050262         3552.686  0.021100330      3527.111 -0.031602951
#> 2:  0.06882734         2927.027  0.013989979      2894.040  0.054837357
#> 3:  0.10570442         3258.389  0.102853985      2748.745  0.002850436
#> 4: -0.06486347         3031.551 -0.062987209      2492.254 -0.001876258
#> 5:  0.08746122         3580.484  0.049451496      3851.801  0.038009725
#> 6: -0.01276267         3915.182  0.009125326      4338.094 -0.021887995
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        3578.260        s1 0.967
#> 2:        2960.015        s2 0.232
#> 3:        3768.034        s3 0.176
#> 4:        3570.847        s4 1.000
#> 5:        3309.168        s5 0.968
#> 6:        3492.270        s6 0.001

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
#>             s1            s2           s3           s4          s5
#> gs1 0.02110033  0.0139899795  0.102853985 -0.062987209  0.04945150
#> gs2 0.03160295 -0.0548373575 -0.002850436  0.001876258 -0.03800973
#> gs3 0.03601390  0.0009200736  0.034431160 -0.041643331  0.01820841
#> gs4 0.04529510  0.1164607591  0.010828368  0.016604712  0.04353669
#> gs5 0.03496842  0.0097605152 -0.012306948 -0.030015126 -0.04707899
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
#>       expr      min       lq     mean   median       uq       max neval
#>  singscore 96.49141 96.50629 97.72522 97.52181 97.70302 100.40354     5
#>   bixverse 18.68955 18.87083 21.46219 20.49903 21.37384  27.87768     5
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
#>     mitch 164.831608 166.881344 171.512672 169.171892 175.135853 181.542661
#>  bixverse   3.580125   3.669587   4.972614   3.809936   5.846095   7.957325
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
