# Pathway activity functions

## Single sample pathway activity functions

An often used algorithm to check pathway activity (or maybe rather
concordant up or down-regulation of genes belonging to a given pathway)
is
[GSVA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7).
There is a corresponding package on
[BioConductor](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)
for this, but `bixverse` offers Rust-accelerated versions of two of the
algorithm.

``` r
library(bixverse)
```

### Data

Similar to the vignette of GSVA, let us create some random data
(Gaussian or Poisson distributed) to mimic different types of
transcriptomics measurements.

``` r
p <- 10000 # number of genes
n <- 100   # number of samples
  
# simulate expression values from a standard Gaussian distribution
X <- matrix(
  rnorm(p * n),
  nrow = p,
  dimnames = list(paste0("g", 1:p), paste0("s", 1:n))
)

# simulate expression values from a Poisson distribution
X1 <- matrix(
  rpois(p * n, lambda = 10),
  nrow = p,
  dimnames = list(paste0("g", 1:p), paste0("s", 1:n))
)
storage.mode(X1) <- "numeric"

# generate some random gene sets
gs <- as.list(sample(10:100, size = 250, replace = TRUE))
# sample gene sets
gs <- lapply(
  gs,
  function(n, p) paste0("g", sample(1:p, size = n, replace = FALSE)),
  p
)
names(gs) <- paste0("gs", 1:length(gs))
```

### GSVA (Gaussian version)

To run the GSVA implementation in bixverse, you can just do

``` r
bixverse_res_gaussian <- calc_gsva(
  exp = X,
  pathways = gs,
  gaussian = TRUE
)

bixverse_res_gaussian[1:5, 1:5]
#>              s1           s2          s3           s4          s5
#> gs1 -0.06895899 -0.071546892  0.05858665 -0.198480721 -0.10222175
#> gs2 -0.09195779  0.038801858  0.02206021  0.091776778 -0.16374353
#> gs3 -0.25206388  0.009658306 -0.03040173 -0.222839353  0.03857107
#> gs4  0.11730874  0.067487905 -0.20165586  0.009264853  0.02935226
#> gs5  0.07854872 -0.123887834 -0.01572874  0.116935674  0.23982384
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1          s2          s3           s4          s5
#> gs1 -0.06895524 -0.07154885  0.05858593 -0.198478818 -0.10221785
#> gs2 -0.09195465  0.03879579  0.02206004  0.091671715 -0.16374657
#> gs3 -0.25204607  0.00965783 -0.03041354 -0.222853452  0.03855260
#> gs4  0.11720804  0.06748753 -0.20165195  0.009253523  0.02946413
#> gs5  0.07855944 -0.12389574 -0.01573449  0.116934616  0.23981520
```

We can appreciated very similar values. There is slight differences in
numerical precision due to optimisations in Rust. However, the
correlations between the two are ≥0.99

``` r
correlations <- diag(
  cor(bixverse_res_gaussian, gsva_res_gaussian)
)

same_res <- all(correlations >= 0.99)

print(
  paste(
    "GSVA and bixverse return (basically) the",
    sprintf("same results for GSVA (Gaussian): %s", same_res)
  )
)
#> [1] "GSVA and bixverse return (basically) the same results for GSVA (Gaussian): TRUE"
```

The `biverse` code is highly optimised, hence, usually yielding faster
results compared to the original `GSVA` code.

``` r
microbenchmark::microbenchmark(
  gsva = {
    gsvaPar <- gsvaParam(X, gs)
    gsva(gsvaPar, verbose = FALSE)
  },
  bixverse = calc_gsva(
    exp = X,
    pathways = gs,
    gaussian = TRUE
  ),
  times = 10L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 2470.9697 2488.8127 2578.6087 2491.4899 2512.4444 3351.1608    10
#>  bixverse  501.4167  501.5525  502.6824  502.7213  503.5358  504.2101    10
```

### GSVA (Poisson version)

To use the version with a Poisson kernel (more appropriate for count
data) is as simple.

``` r
# bixverse expects floats also for the Poisson version
storage.mode(X1) <- "numeric"

bixverse_res_poisson <- calc_gsva(
  exp = X1,
  pathways = gs,
  gaussian = FALSE # uses the Poisson kernel now
)

bixverse_res_poisson[1:5, 1:5]
#>              s1           s2          s3          s4          s5
#> gs1  0.24698656  0.113705234  0.11254477 -0.02482308  0.01390335
#> gs2 -0.01926938  0.219645903  0.24749984 -0.10601103  0.14352004
#> gs3  0.00563511  0.023546120 -0.01558739 -0.09336187  0.44860162
#> gs4 -0.12724825  0.009298187 -0.06946952 -0.06769708 -0.05831338
#> gs5 -0.13849358 -0.151662544 -0.06270434  0.12311371 -0.16182209
```

This is how you would do this in the official GSVA code:

``` r
# update this to use for kcdf = "Poisson"
gsvaParPoisson <- gsvaParam(X1, gs, kcdf = "Poisson")

gsva_res_poisson <- as.matrix(
  gsva(gsvaParPoisson, verbose=FALSE)
)

# Check the correlation between the two
correlations <- diag(
  cor(bixverse_res_poisson, gsva_res_poisson)
)

# basically the same results
same_res <- all(correlations >= 0.99)

print(
  paste(
    "GSVA and bixverse return (basically) the",
    sprintf("same results for GSVA (Poisson): %s", same_res)
  )
)
#> [1] "GSVA and bixverse return (basically) the same results for GSVA (Poisson): TRUE"

gsva_res_poisson[1:5, 1:5]
#>              s1           s2          s3          s4          s5
#> gs1  0.24698656  0.113705234  0.11254477 -0.02482308  0.01390335
#> gs2 -0.01926938  0.219645903  0.24749984 -0.10601103  0.14352004
#> gs3  0.00563511  0.023546120 -0.01558739 -0.09336187  0.44860162
#> gs4 -0.12724825  0.009298187 -0.06946952 -0.06769708 -0.05831338
#> gs5 -0.13849358 -0.151662544 -0.06270434  0.12311371 -0.16182209
```

And speed comparisons for the Poisson kernel version between `GSVA` and
`bixverse`:

``` r
microbenchmark::microbenchmark(
  gsva = {
    gsvaPar <- gsvaParam(X1, gs, kcdf = "Poisson")
    gsva(gsvaPar, verbose = FALSE)
  },
  bixverse = calc_gsva(
    exp = X1,
    pathways = gs,
    gaussian = FALSE
  ),
  times = 10L
)
#> Unit: seconds
#>      expr      min        lq      mean    median        uq       max neval
#>      gsva 18.83928 18.848976 18.866620 18.866146 18.876146 18.912628    10
#>  bixverse  2.71192  2.717477  2.726692  2.720179  2.733643  2.758285    10
```

### single sample GSEA

ssGSEA was first described in [Barbie et
al.](https://www.nature.com/articles/nature08460) and there is an
implementation in `GSVA`. Compared to `GSVA` there is no normalisation
across samples based on the kcdf. `bixverse` also provides a Rust
optimised version of this algorithm.

``` r
bixverse_res_ssgsea <- calc_ssgsea(
  exp = X,
  pathways = gs
)

bixverse_res_ssgsea[1:5, 1:5]
#>              s1          s2          s3         s4          s5
#> gs1  0.07473584  0.07444624 0.189762489 0.01143312  0.02562908
#> gs2  0.03884850  0.12031532 0.154772444 0.13860652 -0.01693133
#> gs3 -0.05141448  0.06852297 0.016029932 0.16402875  0.15117389
#> gs4  0.19583322  0.08284267 0.009791545 0.11028672  0.20345725
#> gs5  0.15452458 -0.01796394 0.091799237 0.19240711  0.25362209
```

This is the way you would run the ssGSEA algorithm in the GSVA package.

``` r
# update this to the ssgsea parameters
ssgseaPar <- ssgseaParam(X, gs)

ssgsea_res <- as.matrix(
  gsva(ssgseaPar, verbose=FALSE)
)

# Check the correlation between the two
correlations <- diag(
  cor(bixverse_res_ssgsea, ssgsea_res)
)

# basically the same results
same_res <- all(correlations >= 0.99)

print(
  paste(
    "GSVA and bixverse return (basically) the",
    sprintf("same results for ssGSEA: %s", same_res)
  )
)
#> [1] "GSVA and bixverse return (basically) the same results for ssGSEA: TRUE"

ssgsea_res[1:5, 1:5]
#>              s1          s2          s3         s4          s5
#> gs1  0.07473584  0.07444624 0.189762489 0.01143312  0.02562908
#> gs2  0.03884850  0.12031532 0.154772444 0.13860652 -0.01693133
#> gs3 -0.05141448  0.06852297 0.016029932 0.16402875  0.15117389
#> gs4  0.19583322  0.08284267 0.009791545 0.11028672  0.20345725
#> gs5  0.15452458 -0.01796394 0.091799237 0.19240711  0.25362209
```

Again, you should observe speed improvements thanks to algorithm
improvements and leveraging Rust’s performance.

``` r
microbenchmark::microbenchmark(
  gsva = {
    ssgseaPar <- ssgseaParam(X, gs)
    gsva(ssgseaPar, verbose = FALSE)
  },
  bixverse = calc_ssgsea(
    exp = X,
    pathways = gs
  ),
  times = 10L
)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>      gsva 633.8024 647.0886 653.6825 650.0772 656.3954 693.9911    10
#>  bixverse 114.9779 115.4705 116.6597 115.7936 116.9732 121.3388    10
```
