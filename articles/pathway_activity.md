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
#>               s1          s2           s3          s4           s5
#> gs1 -0.162827524 -0.14299527  0.063546569  0.19487018 -0.066199459
#> gs2 -0.105560160  0.10791042  0.282454633  0.21947363 -0.005104717
#> gs3 -0.009179827 -0.04213188 -0.046474575  0.02561211 -0.002252421
#> gs4 -0.224960347  0.09139254  0.004257214 -0.10141521  0.084674067
#> gs5  0.171749114 -0.35080509  0.139235375 -0.10358006  0.087328084
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
#>      expr       min        lq      mean    median       uq       max neval
#>      gsva 2588.1481 2591.9126 2600.9195 2595.6772 2607.305 2618.9334     3
#>  bixverse  169.0785  171.0935  172.7515  173.1086  174.588  176.0673     3
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
#>              s1            s2          s3          s4           s5
#> gs1 -0.12615798 -0.0008376684  0.09148418 -0.11199198  0.006449348
#> gs2 -0.10302942 -0.0006631410 -0.10852573 -0.10613093 -0.040500077
#> gs3  0.05504367 -0.0218741958  0.10556842  0.27232615 -0.165617421
#> gs4  0.10735164  0.0294974919 -0.06052016  0.05325169  0.045347909
#> gs5  0.23765793 -0.1662366534 -0.18155616  0.21169901  0.157413039
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
#>      expr        min         lq       mean     median         uq       max
#>      gsva 13022.6878 13062.1662 13087.5302 13101.6446 13119.9515 13138.258
#>  bixverse   260.9036   263.6806   271.8474   266.4575   277.3192   288.181
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
#>               s1          s2         s3          s4         s5
#> gs1  0.032183426  0.04978030 0.16229100  0.14037133 0.03403231
#> gs2  0.038053405  0.14786505 0.23653128  0.13701053 0.11776353
#> gs3  0.053200708  0.14361700 0.08826673  0.14856273 0.06444265
#> gs4 -0.005473408  0.19130939 0.08682039  0.03955608 0.09404664
#> gs5  0.105409186 -0.04645958 0.15125310 -0.02958135 0.21261812
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
#>      gsva 413.46736 418.75244 438.04754 424.03752 450.33763 476.63774     3
#>  bixverse  25.63626  27.77861  29.23459  29.92097  31.03375  32.14653     3
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
#> g1 4033 6636 1618   19 9262
#> g2 4547 6251 2309 2820 3015
#> g3 9412 4728 5074 7333 1987
#> g4 2745 9809 1063 8235 7159
#> g5 6185 2121 4026 1205 6221
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
#>     total_score total_dispersion    up_score up_dispersion   down_score
#>           <num>            <num>       <num>         <num>        <num>
#> 1: -0.006635908         3811.770 -0.04123575      4063.813  0.034599840
#> 2: -0.066003328         3968.185 -0.03480215      3939.274 -0.031201181
#> 3: -0.059507911         3898.873  0.02461603      4564.932 -0.084123940
#> 4:  0.005779659         3041.929  0.04536888      2903.676 -0.039589221
#> 5: -0.038888602         3695.757 -0.02963280      3321.029 -0.009255806
#> 6:  0.020021993         3772.111  0.01817237      3946.687  0.001849625
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        3559.728        s1
#> 2:        3997.096        s2
#> 3:        3232.814        s3
#> 4:        3180.182        s4
#> 5:        4070.484        s5
#> 6:        3597.534        s6
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
#>     total_score total_dispersion    up_score up_dispersion   down_score
#>           <num>            <num>       <num>         <num>        <num>
#> 1: -0.006635908         3811.770 -0.04123575      4063.813  0.034599840
#> 2: -0.066003328         3968.185 -0.03480215      3939.274 -0.031201181
#> 3: -0.059507911         3898.873  0.02461603      4564.932 -0.084123940
#> 4:  0.005779659         3041.929  0.04536888      2903.676 -0.039589221
#> 5: -0.038888602         3695.757 -0.02963280      3321.029 -0.009255806
#> 6:  0.020021993         3772.111  0.01817237      3946.687  0.001849625
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        3559.728        s1 0.976
#> 2:        3997.096        s2 1.000
#> 3:        3232.814        s3 1.000
#> 4:        3180.182        s4 0.039
#> 5:        4070.484        s5 0.972
#> 6:        3597.534        s6 0.007

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
#>              s1          s2           s3          s4           s5
#> gs1 -0.04123575 -0.03480215  0.024616030  0.04536888 -0.029632797
#> gs2 -0.03459984  0.03120118  0.084123940  0.03958922  0.009255806
#> gs3 -0.02625616  0.01212554 -0.004906483  0.02593846 -0.006905822
#> gs4 -0.07690293  0.06012407 -0.005377401 -0.02284033  0.006955029
#> gs5  0.03173336 -0.08128819  0.044207506 -0.06881404  0.065555562
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
#>  singscore 64.80983 64.84016 65.06399 65.05103 65.07987 65.53904     5
#>   bixverse 11.64771 11.90059 12.10648 11.90613 11.96300 13.11496     5
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
#> s1   8.60  0.01 -0.55  0.35
#> s2   7.84  0.40  0.72  2.06
#> s3   8.10 -1.37 -0.54 -1.12
#> s11 -0.42 -1.32  0.22 -0.11
#> s12 -1.72  1.03  0.74  0.19
#> s13  0.65 -0.03  0.96 -0.76
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
#> [1] 0.9992041
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
#>      expr        min         lq       mean     median         uq       max
#>     mitch 129.172759 133.164535 135.723755 134.495400 134.706225 147.07986
#>  bixverse   2.169577   2.319831   3.220868   2.548483   3.999337   5.06711
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
