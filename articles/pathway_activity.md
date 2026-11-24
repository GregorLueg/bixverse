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
#> gs1  0.05635104  0.12046743  0.17078989 -0.10153515  0.08248389
#> gs2  0.02714769  0.02154373  0.04005222 -0.20435932 -0.16290237
#> gs3  0.01035897 -0.14598230 -0.07409544  0.13404508  0.28538150
#> gs4 -0.34082490 -0.20409721 -0.09339093 -0.08006994  0.26680229
#> gs5  0.12936129 -0.31589490 -0.07680595  0.15336488 -0.07906816
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
#>      gsva 2348.3568 2349.8717 2359.0819 2351.3867 2364.4445 2377.5023     3
#>  bixverse  205.5018  206.0758  208.4743  206.6497  209.9605  213.2714     3
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
#> gs1 -0.17655843 -0.02575751  0.07129161  0.07665055  0.06077828
#> gs2  0.05546564  0.08120170 -0.10858561  0.02742134  0.11026263
#> gs3 -0.02556703 -0.05266725 -0.14203124  0.07266064 -0.03646561
#> gs4 -0.31335839 -0.19921767  0.38680010  0.30437723  0.25911253
#> gs5  0.05949218  0.05147870 -0.28419250 -0.37843784 -0.12616375
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
#>      expr        min         lq       mean     median         uq        max
#>      gsva 12493.7064 12502.8787 12506.9360 12512.0509 12513.5508 12515.0506
#>  bixverse   299.6463   299.9096   300.4337   300.1729   300.8274   301.4819
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
#>             s1          s2           s3         s4         s5
#> gs1  0.1483540  0.18818806  0.243765565 0.10211428 0.18749134
#> gs2  0.1155501  0.15804823  0.160009567 0.04582334 0.03334051
#> gs3  0.1151616  0.05328887  0.041567412 0.17909917 0.24178154
#> gs4 -0.0737504  0.01994511 -0.007312245 0.19690389 0.25686273
#> gs5  0.2648578 -0.10301076  0.063050263 0.21927844 0.08259099
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
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 503.67924 504.09681 508.07552 504.51437 510.27367 516.03296     3
#>  bixverse  34.12563  36.16291  37.08364  38.20019  38.56264  38.92508     3
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
#> g1 6006 2523 6144 6930 4524
#> g2 7350 6436 3173 3141 5487
#> g3 8032  633 3304 3248 7002
#> g4 8855 7322 3094  399 8832
#> g5 3828 1149   42 4681 9357
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
#> 1:  0.02172272         3567.141  0.016962626      3282.481  0.004760097
#> 2:  0.01381126         3832.156  0.027905051      4160.182 -0.014093793
#> 3:  0.03583345         3473.737  0.058096970      3889.607 -0.022263520
#> 4:  0.02342985         3779.524 -0.017079798      3762.103  0.040509645
#> 5:  0.05528355         3445.938  0.025362626      3696.127  0.029920925
#> 6: -0.02635067         3662.769 -0.001177778      3838.457 -0.025172895
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        3851.801        s1
#> 2:        3504.130        s2
#> 3:        3057.867        s3
#> 4:        3796.944        s4
#> 5:        3195.749        s5
#> 6:        3487.080        s6
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
#> 1:  0.02172272         3567.141  0.016962626      3282.481  0.004760097
#> 2:  0.01381126         3832.156  0.027905051      4160.182 -0.014093793
#> 3:  0.03583345         3473.737  0.058096970      3889.607 -0.022263520
#> 4:  0.02342985         3779.524 -0.017079798      3762.103  0.040509645
#> 5:  0.05528355         3445.938  0.025362626      3696.127  0.029920925
#> 6: -0.02635067         3662.769 -0.001177778      3838.457 -0.025172895
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        3851.801        s1 0.085
#> 2:        3504.130        s2 0.001
#> 3:        3057.867        s3 0.001
#> 4:        3796.944        s4 0.711
#> 5:        3195.749        s5 0.203
#> 6:        3487.080        s6 1.000

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
#>               s1          s2          s3          s4          s5
#> gs1  0.016962626  0.02790505  0.05809697 -0.01707980  0.02536263
#> gs2 -0.004760097  0.01409379  0.02226352 -0.04050965 -0.02992092
#> gs3  0.002709073 -0.04632422 -0.04409180  0.03524588  0.07409828
#> gs4 -0.081800235 -0.04569969 -0.06539155  0.01896226  0.09664960
#> gs5  0.077984962 -0.09536040 -0.01659549  0.05410526 -0.03918797
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
#>  singscore 64.49255 64.77806 65.32824 65.72029 65.73363 65.91669     5
#>   bixverse 13.42732 15.63357 15.41740 15.77794 16.03292 16.21526     5
```

## ULM and MLM

The rank-based scores above treat a gene set as a bag of genes.
Transcription factor regulons (think
[CollecTRI](https://academic.oup.com/nar/article/52/D1/D1029/7386439) or
[DoRothEA](https://genome.cshlp.org/content/29/8/1363)) carry a sign as
well: some targets get activated, others repressed.
[decoupler](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613)
handles that with a linear model per sample: regress the expression of
all genes on the weights of a source’s targets, and take the t-value of
the slope as the activity. The univariate version (ULM) fits one model
per source; the multivariate one (MLM) puts all sources into one model,
so overlapping regulons compete for the signal.
[`calc_lm_activity()`](https://gregorlueg.github.io/bixverse/reference/calc_lm_activity.md)
is the Rust version of both.

The network comes in decoupler’s long format: `source`, `target` and an
optional `weight`. We turn the random gene sets into signed regulons and
plant an activation of the first one in the first ten samples.

``` r

net <- data.table::rbindlist(
  purrr::imap(gs[1:50], \(targets, source) {
    data.table::data.table(
      source = source,
      target = targets,
      weight = sample(c(-1, 1), length(targets), replace = TRUE)
    )
  })
)

X_lm <- X
gs1_net <- net[source == "gs1"]
X_lm[gs1_net$target, 1:10] <- X_lm[gs1_net$target, 1:10] + gs1_net$weight
```

``` r

ulm_res <- calc_lm_activity(
  mat = X_lm,
  net = net,
  lm_params = params_lm_activity(method = "ulm", p_values = TRUE)
)

# planted samples first, the rest after
round(ulm_res$es[c(1:3, 11:13), 1:4], 2)
#>       gs1  gs10  gs11  gs12
#> s1  12.07  0.66  1.20 -0.28
#> s2   8.91  0.53  2.41  1.68
#> s3   9.54  1.37 -0.03 -0.58
#> s11  1.94 -1.55 -0.80  0.29
#> s12  2.04  1.75 -0.89  1.63
#> s13  0.75  0.15 -0.33 -0.32
```

`es` holds the t-values (`stat = "coef"` gives you the slope instead),
`pvals` the two-sided p-values and `fdr` the Benjamini-Hochberg
adjustment per sample across sources, which is what decoupler does.
Sources with fewer than `tmin` targets in the data (default 5) are
dropped, same as decoupler.

``` r

mlm_res <- calc_lm_activity(
  mat = X_lm,
  net = net,
  lm_params = params_lm_activity(method = "mlm", p_values = TRUE)
)

cor(ulm_res$es[, "gs1"], mlm_res$es[, "gs1"])
#> [1] 0.999755
```

A few things to be aware of:

- MLM throws an error on collinear sources, i.e., duplicated or
  proportional regulons. Merge or drop them first.
- Rows with `NA` are dropped (with a warning) before scoring.
- The p-values treat genes as independent. They are not, so read them as
  a ranking aid rather than gospel.
- A named vector (say the t statistics of a contrast) works as well and
  gives you a `data.table` back.

### Comparison with decoupleR

``` r

dc_ulm <- decoupleR::run_ulm(
  mat = X_lm,
  network = net,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 5L
)

dc_ulm_wide <- data.table::dcast(
  data.table::as.data.table(dc_ulm),
  condition ~ source,
  value.var = "score"
)
dc_ulm_mat <- as.matrix(dc_ulm_wide[, -1])
rownames(dc_ulm_mat) <- dc_ulm_wide$condition

max(abs(dc_ulm_mat[rownames(ulm_res$es), colnames(ulm_res$es)] - ulm_res$es))
```

And the speed difference:

``` r

microbenchmark::microbenchmark(
  decoupleR_ulm = decoupleR::run_ulm(
    mat = X_lm,
    network = net,
    .source = "source",
    .target = "target",
    .mor = "weight"
  ),
  bixverse_ulm = calc_lm_activity(mat = X_lm, net = net),
  decoupleR_mlm = decoupleR::run_mlm(
    mat = X_lm,
    network = net,
    .source = "source",
    .target = "target",
    .mor = "weight"
  ),
  bixverse_mlm = calc_lm_activity(
    mat = X_lm,
    net = net,
    lm_params = params_lm_activity(method = "mlm")
  ),
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
#>      expr        min        lq       mean   median         uq        max neval
#>     mitch 112.365064 117.86223 120.941509 118.1219 123.614481 132.743914     5
#>  bixverse   2.727863   2.75401   3.739062   3.0787   4.692355   5.442384     5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
