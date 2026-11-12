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
#>               s1          s2         s3          s4          s5
#> gs1  0.009293610 -0.29871006 0.06977003 -0.04612870  0.14043544
#> gs2 -0.028602657  0.16253117 0.06744198 -0.03472199 -0.09320403
#> gs3  0.103004671  0.09122989 0.04839992 -0.05013398  0.20757637
#> gs4 -0.007525981 -0.02252031 0.11368999  0.10205669 -0.15187084
#> gs5 -0.054715083 -0.07994525 0.20715640 -0.13225299 -0.06335602
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
#>      gsva 2958.7746 2967.0418 2979.4814 2975.3090 2989.8347 3004.3604     3
#>  bixverse  225.2219  226.9961  227.5883  228.7703  228.7716  228.7728     3
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
#>              s1          s2           s3          s4          s5
#> gs1 -0.15486794  0.10811213  0.042165974 -0.27444985 -0.08189733
#> gs2 -0.04630420  0.09568645  0.198676732 -0.17817616  0.07110626
#> gs3 -0.19122926  0.07873679 -0.098723565  0.13148680  0.01224628
#> gs4  0.10444568 -0.04123473  0.005706957  0.02007722  0.02453384
#> gs5 -0.07699979 -0.07954422 -0.015529801 -0.07392899 -0.00855247
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
#>      gsva 15512.5853 15516.1839 15529.9031 15519.7824 15538.5620 15557.3415
#>  bixverse   357.8424   359.8259   365.3963   361.8094   369.1733   376.5372
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
#> gs1 0.08707479 0.03272963 0.11936042 0.11105499 0.21254804
#> gs2 0.12534760 0.15961035 0.08223456 0.09435023 0.05631715
#> gs3 0.12397934 0.13119054 0.14412321 0.07694772 0.21477382
#> gs4 0.07647943 0.11206578 0.13217191 0.11372099 0.01847546
#> gs5 0.05812219 0.07718779 0.20429049 0.01871027 0.06849295
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
#>      expr       min        lq      mean    median       uq       max neval
#>      gsva 533.49704 535.77549 539.95073 538.05394 543.1776 548.30121     3
#>  bixverse  33.38343  34.25034  35.77408  35.11725  36.9694  38.82154     3
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
#> g1 1544 3407 6514 5888 4498
#> g2 5212 8876 2423 7751 5539
#> g3 9253 3465 4327 9703 2367
#> g4 1159 3100 7667 6632  681
#> g5 2611 2998 1199 7119 4849
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
#>     total_score total_dispersion     up_score up_dispersion    down_score
#>           <num>            <num>        <num>         <num>         <num>
#> 1:  0.000226477         4196.876  0.001821037      4249.879 -0.0015945602
#> 2: -0.107101761         4015.999 -0.068134356      4656.112 -0.0389674054
#> 3:  0.012434809         3702.799  0.012120756      4103.843  0.0003140524
#> 4:  0.005037007         3939.274  0.002923684      4200.212  0.0021133237
#> 5:  0.089757107         3952.247  0.062491647      4616.082  0.0272654599
#> 6:  0.044177166         3866.256  0.025669941      4964.494  0.0185072244
#>    down_dispersion sample_id
#>              <num>    <char>
#> 1:        4143.873        s1
#> 2:        3375.885        s2
#> 3:        3301.755        s3
#> 4:        3678.336        s4
#> 5:        3288.412        s5
#> 6:        2768.018        s6
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
#>     total_score total_dispersion     up_score up_dispersion    down_score
#>           <num>            <num>        <num>         <num>         <num>
#> 1:  0.000226477         4196.876  0.001821037      4249.879 -0.0015945602
#> 2: -0.107101761         4015.999 -0.068134356      4656.112 -0.0389674054
#> 3:  0.012434809         3702.799  0.012120756      4103.843  0.0003140524
#> 4:  0.005037007         3939.274  0.002923684      4200.212  0.0021133237
#> 5:  0.089757107         3952.247  0.062491647      4616.082  0.0272654599
#> 6:  0.044177166         3866.256  0.025669941      4964.494  0.0185072244
#>    down_dispersion sample_id  pval
#>              <num>    <char> <num>
#> 1:        4143.873        s1 0.001
#> 2:        3375.885        s2 1.000
#> 3:        3301.755        s3 1.000
#> 4:        3678.336        s4 0.033
#> 5:        3288.412        s5 0.001
#> 6:        2768.018        s6 0.933

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
#>               s1           s2            s3           s4          s5
#> gs1  0.001821037 -0.068134356  0.0121207565  0.002923684  0.06249165
#> gs2  0.001594560  0.038967405 -0.0003140524 -0.002113324 -0.02726546
#> gs3  0.008025842  0.022423952  0.0217300676 -0.012477337  0.06680460
#> gs4 -0.012599705  0.005730393  0.0236985299  0.018608530 -0.05142769
#> gs5 -0.015363520 -0.017619318  0.0571377150 -0.045155603 -0.02064383
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
#>  singscore 76.59327 77.55136 77.52306 77.77990 77.81225 77.87850     5
#>   bixverse 14.15412 14.36551 14.37317 14.41455 14.46302 14.46866     5
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
#> s1   4.71 -0.74 -0.74 -0.98
#> s2   6.70 -0.07 -0.68 -1.07
#> s3   3.65  0.76  0.08  1.40
#> s11 -0.08 -0.20  0.08 -0.61
#> s12 -1.85  0.69  0.00 -0.21
#> s13  0.18 -0.02  0.69 -1.19
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
#> [1] 0.9997787
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
#>     mitch 147.918207 150.000301 154.651165 150.796601 158.410252 166.130462
#>  bixverse   2.664316   2.681061   3.605618   2.799987   4.187682   5.695045
#>  neval
#>      5
#>      5
```

As with the GSVA and ssGSEA implementations, you should observe
meaningful speed improvements from Rust, particularly as the number of
gene sets or contrasts grows and the more cores/oomph your system has.
