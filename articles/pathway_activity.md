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
#> gs1  0.070699716  0.24630664 -0.061413574 -0.11923662 -0.162808931
#> gs2 -0.060057565 -0.04325684 -0.063622669  0.09660418 -0.059086540
#> gs3  0.007450875 -0.03024081  0.111842308 -0.11288680  0.049043434
#> gs4 -0.147745032  0.20495200 -0.001324612 -0.08988726  0.029333428
#> gs5  0.015130146  0.09196079 -0.133439495  0.04130469 -0.007387446
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
#>      expr      min        lq      mean    median        uq      max neval
#>      gsva 2952.676 2955.0143 2962.1652 2957.3525 2966.9098 2976.467     3
#>  bixverse  225.291  226.3476  226.8064  227.4043  227.5642  227.724     3
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
#> gs1 -0.00118558 -0.04826283  0.08271244 -0.008508838  0.02148845
#> gs2  0.19879189  0.12793483 -0.12120967 -0.195681004  0.04726362
#> gs3 -0.10323993 -0.02025866  0.10515224  0.083007988  0.17278013
#> gs4  0.06414662  0.04450705  0.06484290  0.059778216 -0.14512430
#> gs5 -0.11058538 -0.05783036  0.00419940 -0.049852569  0.03254594
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
#>      gsva 15466.6979 15470.7975 15476.7913 15474.8972 15481.8381 15488.7790
#>  bixverse   344.1202   345.7618   352.2213   347.4034   356.2719   365.1404
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
#>             s1         s2         s3         s4         s5
#> gs1 0.10396213 0.22974634 0.05447958 0.10871633 0.04833756
#> gs2 0.03043592 0.10767493 0.10304002 0.18185985 0.05853864
#> gs3 0.10118229 0.07397841 0.12400023 0.06478454 0.16938805
#> gs4 0.05956131 0.22237959 0.06720025 0.06150574 0.09182176
#> gs5 0.09573253 0.15485268 0.04466422 0.11203095 0.09815160
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
#>      gsva 533.97955 535.20619 541.88892 536.43284 545.84360 555.25437     3
#>  bixverse  37.67757  38.27871  43.55546  38.87986  46.49441  54.10897     3
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
#> g1 2709 7862 5802 3999 1205
#> g2 1056 5766 5678 6352 3810
#> g3 2418  791 4200 1629 1777
#> g4 5233 8449 3249 7993   32
#> g5 6279 4236 6422 4844 1967
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
#>     total_score total_dispersion     up_score up_dispersion   down_score
#>           <num>            <num>        <num>         <num>        <num>
#> 1:  0.029780005         3556.021 -0.006354380      3917.035  0.036134385
#> 2:  0.072532842         3216.876  0.077659401      3192.043 -0.005126560
#> 3: -0.029179939         3730.598 -0.033315783      3484.115  0.004135844
#> 4: -0.055133993         3676.854 -0.008466571      4108.291 -0.046667421
#> 5: -0.022951815         3753.949 -0.045623476      4194.282  0.022671661
#> 6:  0.006263995         3703.170 -0.007577717      4032.678  0.013841712
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        3195.008        s1
#> 2:        3241.710        s2
#> 3:        3977.080        s3
#> 4:        3245.416        s4
#> 5:        3313.616        s5
#> 6:        3373.661        s6
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
#>     total_score total_dispersion     up_score up_dispersion   down_score
#>           <num>            <num>        <num>         <num>        <num>
#> 1:  0.029780005         3556.021 -0.006354380      3917.035  0.036134385
#> 2:  0.072532842         3216.876  0.077659401      3192.043 -0.005126560
#> 3: -0.029179939         3730.598 -0.033315783      3484.115  0.004135844
#> 4: -0.055133993         3676.854 -0.008466571      4108.291 -0.046667421
#> 5: -0.022951815         3753.949 -0.045623476      4194.282  0.022671661
#> 6:  0.006263995         3703.170 -0.007577717      4032.678  0.013841712
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        3195.008        s1 0.001
#> 2:        3241.710        s2 0.001
#> 3:        3977.080        s3 1.000
#> 4:        3245.416        s4 1.000
#> 5:        3313.616        s5 0.892
#> 6:        3373.661        s6 0.001

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
#>                s1          s2           s3           s4            s5
#> gs1 -0.0063543804  0.07765940 -0.033315783 -0.008466571 -0.0456234757
#> gs2 -0.0361343851  0.00512656 -0.004135844  0.046667421 -0.0226716611
#> gs3  0.0001072799 -0.01328483  0.008089114 -0.022796976  0.0301424929
#> gs4 -0.0277948712  0.06749541 -0.009659317 -0.023272214 -0.0174807086
#> gs5  0.0019324202  0.01911183 -0.037118548  0.006896273  0.0009509614
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
#>       expr      min      lq     mean   median       uq      max neval
#>  singscore 76.52818 76.6906 77.06302 76.97281 77.49955 77.62396     5
#>   bixverse 13.58382 14.0240 14.83336 14.31970 14.98647 17.25281     5
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
#> s1   7.88 -0.52 -1.27  2.02
#> s2   8.01 -0.37 -0.46 -0.99
#> s3   9.55 -0.25  1.02 -0.31
#> s11 -0.25 -2.34  0.55  0.02
#> s12 -0.57 -1.23  0.06  1.54
#> s13  1.43 -1.06  0.62 -0.86
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
#> [1] 0.9996212
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
#>      expr        min       lq      mean     median       uq        max neval
#>     mitch 144.703150 146.5434 154.81524 151.701162 160.1683 170.960147     5
#>  bixverse   2.707881   2.7208   3.86193   2.752938   4.4480   6.680032     5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
