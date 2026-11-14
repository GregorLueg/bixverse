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
#> gs1  0.151660790  0.15090522 -0.21954131 -0.13206454 -0.11402024
#> gs2 -0.312775008  0.26470172  0.19100363  0.10107066  0.08371026
#> gs3  0.344598187  0.53676368 -0.17103953 -0.07744563 -0.04201189
#> gs4  0.141634180  0.22107595  0.05313029 -0.12979672 -0.12457845
#> gs5 -0.003180515 -0.03735953  0.03455718 -0.18735479  0.13472310
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
#>      gsva 3391.2638 3398.2068 3411.6994 3405.150 3421.9172 3438.6845     3
#>  bixverse  262.6093  263.2967  267.8489  263.984  270.4688  276.9535     3
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
#>             s1           s2          s3          s4           s5
#> gs1 -0.1832271  0.093012157  0.10477476  0.13669297 -0.124568171
#> gs2  0.2911998  0.007746107  0.08234871 -0.12439929  0.107213116
#> gs3 -0.3197658 -0.026634298 -0.27708862 -0.06036900 -0.209775914
#> gs4 -0.3341146  0.125733016  0.01413535  0.06486563  0.004427635
#> gs5 -0.1158344 -0.021105012  0.03737275 -0.17065933  0.160874536
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
#>      gsva 19479.082 19491.3996 19499.6215 19503.7170 19509.8913 19516.0655
#>  bixverse   413.905   417.3874   421.9568   420.8698   425.9827   431.0956
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
#>              s1        s2          s3           s4         s5
#> gs1  0.16814269 0.1567567 -0.01763440 -0.015419819 0.01336638
#> gs2 -0.02687646 0.1163365  0.15271729  0.177148892 0.13036791
#> gs3  0.25559946 0.4917374 -0.04050507  0.020693620 0.13436826
#> gs4  0.25050249 0.2295795  0.18279925  0.004007419 0.10061898
#> gs5  0.11727349 0.1380222  0.09565550  0.089161165 0.13820305
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
#>      expr       min        lq      mean   median         uq        max neval
#>      gsva 649.93431 650.82548 940.30172 651.7167 1085.48542 1519.25420     3
#>  bixverse  38.69379  39.55609  40.91291  40.4184   42.02247   43.62654     3
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
#> g1 8576 5856 9513  298 7288
#> g2 1257 5468 6806 8249 9791
#> g3 4616 8666 6543 9945 2043
#> g4 3131 8480 4762 8402 4306
#> g5 5910 9515 8158 9134 8707
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
#>    total_score total_dispersion      up_score up_dispersion  down_score
#>          <num>            <num>         <num>         <num>       <num>
#> 1:  0.11300371         3310.280  0.0424600615      3510.061  0.07054365
#> 2:  0.01665425         3306.574  0.0326497132      3352.905 -0.01599546
#> 3: -0.11333160         3725.409 -0.0708384492      3624.221 -0.04249315
#> 4: -0.08164499         4213.926 -0.0677533140      3378.850 -0.01389167
#> 5: -0.07478188         3997.837 -0.0607176798      3758.397 -0.01406420
#> 6: -0.05562065         3894.425  0.0001352673      4360.333 -0.05575592
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        3110.499        s1
#> 2:        3260.242        s2
#> 3:        3826.596        s3
#> 4:        5049.002        s4
#> 5:        4237.277        s5
#> 6:        3428.518        s6
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
#>    total_score total_dispersion      up_score up_dispersion  down_score
#>          <num>            <num>         <num>         <num>       <num>
#> 1:  0.11300371         3310.280  0.0424600615      3510.061  0.07054365
#> 2:  0.01665425         3306.574  0.0326497132      3352.905 -0.01599546
#> 3: -0.11333160         3725.409 -0.0708384492      3624.221 -0.04249315
#> 4: -0.08164499         4213.926 -0.0677533140      3378.850 -0.01389167
#> 5: -0.07478188         3997.837 -0.0607176798      3758.397 -0.01406420
#> 6: -0.05562065         3894.425  0.0001352673      4360.333 -0.05575592
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        3110.499        s1 0.016
#> 2:        3260.242        s2 1.000
#> 3:        3826.596        s3 1.000
#> 4:        5049.002        s4 1.000
#> 5:        4237.277        s5 1.000
#> 6:        3428.518        s6 0.993

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
#>               s1         s2           s3           s4          s5
#> gs1  0.042460061 0.03264971 -0.070838449 -0.067753314 -0.06071768
#> gs2 -0.070543645 0.01599546  0.042493154  0.013891672  0.01406420
#> gs3  0.100894696 0.23364128 -0.061389603 -0.070404891 -0.01455019
#> gs4  0.073807284 0.07298977  0.034717188 -0.055760362 -0.01357375
#> gs5  0.008178579 0.01229456 -0.009097557 -0.008443918  0.02876498
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
#>       expr      min       lq      mean   median       uq       max neval
#>  singscore 88.24352 89.71392  89.91574 90.13021 90.55346  90.93757     5
#>   bixverse 14.35271 17.11310 201.35390 17.31225 18.69114 939.30030     5
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
#>      gs1  gs10  gs11  gs12
#> s1  8.41  0.56 -0.11  0.76
#> s2  8.39 -0.79  0.53 -1.78
#> s3  7.77  0.21  0.96  1.58
#> s11 0.37  0.01 -0.18 -1.03
#> s12 0.50 -0.44 -0.67  0.95
#> s13 0.41  0.08  0.05  0.36
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
#> [1] 0.9997414
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
#>      expr        min         lq       mean     median         uq        max
#>     mitch 196.109617 196.564656 199.054868 196.619338 201.632982 204.347747
#>  bixverse   3.542958   3.714681   4.798825   4.008469   5.233173   7.494843
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
