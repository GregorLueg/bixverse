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
#>              s1           s2          s3          s4           s5
#> gs1 -0.01544329  0.008398733  0.09405273 -0.02405362 -0.046604377
#> gs2 -0.28626351 -0.193715431  0.32416034  0.13540815 -0.003743303
#> gs3  0.11184395 -0.071052103 -0.09339799 -0.18353478 -0.139473146
#> gs4 -0.09882816 -0.004371325  0.05003732  0.01575666 -0.068463801
#> gs5 -0.08510240  0.083757712  0.14028398 -0.25462415 -0.122935178
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
#>      expr       min        lq      mean    median        uq      max neval
#>      gsva 3161.6541 3187.2178 3199.0115 3212.7816 3217.6901 3222.599     3
#>  bixverse  236.1802  237.3205  239.1197  238.4608  240.5894  242.718     3
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
#> gs1  0.05645088  0.01539904  0.01759720 -0.04009022 -0.01053430
#> gs2  0.13299204 -0.22463535 -0.26928040 -0.10809041 -0.03323701
#> gs3 -0.07951393  0.21022153 -0.24091346 -0.17222717 -0.05191636
#> gs4 -0.01315572  0.06452488 -0.06160873  0.21191037  0.02987137
#> gs5 -0.11875192 -0.37465047 -0.10844028  0.14202819 -0.32907692
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
#> Unit: milliseconds
#>      expr       min         lq       mean     median         uq        max
#>      gsva 16483.574 16496.9463 16510.6173 16510.3183 16524.1389 16537.9594
#>  bixverse   367.762   371.0031   373.1649   374.2441   375.8663   377.4885
#>  neval
#>      3
#>      3
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
#>              s1          s2         s3          s4         s5
#> gs1  0.07071095  0.14225932 0.15678140  0.07686019 0.07749766
#> gs2 -0.07676296 -0.02323372 0.27735071  0.26074213 0.09755611
#> gs3  0.07829024  0.14755386 0.04138406  0.07820624 0.03122962
#> gs4  0.09643127  0.08595270 0.11009523  0.11876077 0.07459953
#> gs5  0.09044204  0.13552273 0.15579406 -0.03488796 0.13428451
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
#>      expr       min       lq      mean    median        uq       max neval
#>      gsva 602.31910 610.7995 614.46143 619.27982 620.53260 621.78538     3
#>  bixverse  55.49415  56.2547  56.82294  57.01525  57.48733  57.95942     3
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
#> g1 8481 3700 4491 5807 2468
#> g2 2625 4541 2205 2447 3966
#> g3 8525 9710 3491  243 4887
#> g4 3874 3636 5610 9458 7739
#> g5 6716 9573  205 1584 2841
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
#>    total_score total_dispersion     up_score up_dispersion  down_score
#>          <num>            <num>        <num>         <num>       <num>
#> 1:  0.08649271         3373.661 -0.028653285      3885.900  0.11514599
#> 2:  0.07550231         3900.726  0.007738554      4146.838  0.06776376
#> 3: -0.07667861         3880.711  0.018295088      4047.504 -0.09497370
#> 4: -0.07834642         3137.928 -0.011672032      3189.077 -0.06667439
#> 5:  0.01665420         3653.132 -0.015616975      3623.480  0.03227117
#> 6: -0.01363445         3693.162 -0.027815229      3598.276  0.01418078
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        2861.422        s1
#> 2:        3654.614        s2
#> 3:        3713.919        s3
#> 4:        3086.778        s4
#> 5:        3682.784        s5
#> 6:        3788.049        s6
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
#>    total_score total_dispersion     up_score up_dispersion  down_score
#>          <num>            <num>        <num>         <num>       <num>
#> 1:  0.08649271         3373.661 -0.028653285      3885.900  0.11514599
#> 2:  0.07550231         3900.726  0.007738554      4146.838  0.06776376
#> 3: -0.07667861         3880.711  0.018295088      4047.504 -0.09497370
#> 4: -0.07834642         3137.928 -0.011672032      3189.077 -0.06667439
#> 5:  0.01665420         3653.132 -0.015616975      3623.480  0.03227117
#> 6: -0.01363445         3693.162 -0.027815229      3598.276  0.01418078
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        2861.422        s1 0.019
#> 2:        3654.614        s2 0.001
#> 3:        3713.919        s3 0.972
#> 4:        3086.778        s4 1.000
#> 5:        3682.784        s5 0.001
#> 6:        3788.049        s6 0.173

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
#>               s1           s2          s3           s4           s5
#> gs1 -0.028653285  0.007738554  0.01829509 -0.011672032 -0.015616975
#> gs2 -0.115145992 -0.067763756  0.09497370  0.066674390 -0.032271173
#> gs3  0.000271572  0.010668530 -0.02646604 -0.024342371 -0.046616419
#> gs4 -0.013467723 -0.013918148  0.01192000  0.009335747 -0.021907503
#> gs5  0.001561624  0.026503040  0.03734713 -0.065430366 -0.002605182
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
#>  singscore 77.38631 77.69995 78.12279 78.34926 78.46742 78.71099     5
#>   bixverse 12.83084 13.10351 13.35220 13.22978 13.63866 13.95821     5
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
#>     mitch 153.663206 154.125381 158.402364 155.471845 161.413185 167.338201
#>  bixverse   3.256484   3.337464   4.313311   3.692579   4.906007   6.374021
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
