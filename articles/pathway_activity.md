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
#>              s1          s2          s3          s4         s5
#> gs1 -0.15928645  0.26795373  0.04438474  0.09671402 0.02991247
#> gs2 -0.00859009 -0.03301556 -0.12917105 -0.14630861 0.04725189
#> gs3 -0.09183870 -0.22465920  0.23222821 -0.12199500 0.14470356
#> gs4 -0.01787853  0.12569636 -0.07943017 -0.10429838 0.10460935
#> gs5 -0.23644842  0.07704799  0.12638993  0.04197266 0.22785429
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>               s1          s2          s3          s4         s5
#> gs1 -0.159275978  0.26805369  0.04438032  0.09672382 0.02990255
#> gs2 -0.008590587 -0.03302900 -0.12906046 -0.14631606 0.04726060
#> gs3 -0.091834134 -0.22465455  0.23223756 -0.12199500 0.14472519
#> gs4 -0.017878006  0.12566689 -0.07941964 -0.10428968 0.10472379
#> gs5 -0.236478699  0.07705751  0.12640330  0.04196244 0.22784787
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
#>      gsva 2488.7947 2492.4342 2693.8834 2508.6248 3134.3301 3139.7895    10
#>  bixverse  580.1808  580.1863  583.4039  582.0213  585.1674  593.4229    10
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
#>              s1          s2          s3          s4          s5
#> gs1  0.07284704  0.03662440  0.11440266  0.01013988  0.02558805
#> gs2 -0.01162582 -0.08673113 -0.02672428  0.03251351 -0.01902238
#> gs3  0.11132842  0.12414347  0.16951170 -0.19409656 -0.36549243
#> gs4 -0.12977008 -0.01305005  0.22092840 -0.03202707 -0.05511910
#> gs5  0.21342239 -0.02838233  0.18356675  0.29247711 -0.15647927
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
#>              s1          s2          s3          s4          s5
#> gs1  0.07284704  0.03662440  0.11440266  0.01013988  0.02558805
#> gs2 -0.01162582 -0.08673113 -0.02672428  0.03251351 -0.01902238
#> gs3  0.11132842  0.12414347  0.16951170 -0.19409656 -0.36549243
#> gs4 -0.12977008 -0.01305005  0.22092840 -0.03202707 -0.05511910
#> gs5  0.21342239 -0.02838233  0.18356675  0.29247711 -0.15647927
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
#>      expr       min        lq     mean    median        uq       max neval
#>      gsva 18.854600 18.859782 18.93452 18.869425 18.886151 19.507765    10
#>  bixverse  2.764678  2.768701  2.77370  2.770128  2.772898  2.808917    10
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
#>               s1         s2         s3         s4         s5
#> gs1  0.063090032 0.24929449 0.15599893 0.11356530 0.09390885
#> gs2  0.132050576 0.09493609 0.06631275 0.03166304 0.10893233
#> gs3 -0.005542842 0.07330195 0.32420435 0.03788568 0.16882813
#> gs4  0.107914830 0.19409809 0.02824002 0.02681649 0.12644325
#> gs5  0.090007140 0.07180718 0.15003539 0.14899284 0.18404771
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
#>               s1         s2         s3         s4         s5
#> gs1  0.063090032 0.24929449 0.15599893 0.11356530 0.09390885
#> gs2  0.132050576 0.09493609 0.06631275 0.03166304 0.10893233
#> gs3 -0.005542842 0.07330195 0.32420435 0.03788568 0.16882813
#> gs4  0.107914830 0.19409809 0.02824002 0.02681649 0.12644325
#> gs5  0.090007140 0.07180718 0.15003539 0.14899284 0.18404771
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
#>      expr       min        lq      mean    median        uq      max neval
#>      gsva 602.69528 606.80449 615.81146 617.42019 620.89029 634.4618    10
#>  bixverse  97.41221  98.33761  98.80699  98.64975  99.31794 100.6047    10
```
