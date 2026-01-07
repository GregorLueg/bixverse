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
#>             s1          s2          s3          s4          s5
#> gs1 0.02531758  0.08710263 -0.12670145  0.16703882  0.06078357
#> gs2 0.16386204 -0.15163228  0.14733795 -0.02592786 -0.20854276
#> gs3 0.03786542  0.10412821  0.09786149  0.05780583 -0.17334996
#> gs4 0.05105737  0.04678924 -0.03830777  0.13725562  0.18453084
#> gs5 0.01372749 -0.17837961 -0.08835323  0.13115592  0.12604453
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>             s1          s2          s3          s4          s5
#> gs1 0.02531758  0.08711483 -0.12669044  0.16704674  0.06078357
#> gs2 0.16384238 -0.15152069  0.14745230 -0.02592334 -0.20853934
#> gs3 0.03786670  0.10412841  0.09786407  0.05779633 -0.17335915
#> gs4 0.05107059  0.04667972 -0.03829334  0.13722062  0.18451889
#> gs5 0.01361604 -0.17838176 -0.08835900  0.13105661  0.12604483
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
#>      gsva 2497.3744 2503.8112 2694.6036 2513.0099 2556.5217 3434.6298    10
#>  bixverse  473.0985  475.6385  477.0952  477.1221  477.9504  481.7166    10
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
#>               s1         s2          s3         s4          s5
#> gs1 -0.119089944  0.0610531  0.19802076  0.0508559 -0.24404625
#> gs2  0.024612460  0.1355144 -0.08305054 -0.1881348 -0.11231032
#> gs3  0.052576086  0.0165992  0.08925749 -0.1719681  0.09611489
#> gs4  0.006948947  0.1511685 -0.07948645 -0.1307248  0.08987483
#> gs5  0.035953880 -0.1186225 -0.04051918  0.1138275  0.02144087
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
#>               s1         s2          s3         s4          s5
#> gs1 -0.119089944  0.0610531  0.19802076  0.0508559 -0.24404625
#> gs2  0.024612460  0.1355144 -0.08305054 -0.1881348 -0.11231032
#> gs3  0.052576086  0.0165992  0.08925749 -0.1719681  0.09611489
#> gs4  0.006948947  0.1511685 -0.07948645 -0.1307248  0.08987483
#> gs5  0.035953880 -0.1186225 -0.04051918  0.1138275  0.02144087
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
#>      expr       min        lq      mean   median        uq       max neval
#>      gsva 18.862805 18.876114 19.003102 18.92172 18.972255 19.809561    10
#>  bixverse  2.683554  2.684942  2.693027  2.68893  2.693078  2.734557    10
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
#>            s1           s2         s3         s4          s5
#> gs1 0.1322743  0.122216769 0.06321831 0.15099473  0.18020272
#> gs2 0.1878729 -0.002001961 0.22666278 0.08462068 -0.08660604
#> gs3 0.1375218  0.131517876 0.13547157 0.11075140  0.08107227
#> gs4 0.1229140  0.169589360 0.10792946 0.14134695  0.23096335
#> gs5 0.1431006  0.009082125 0.09946726 0.21310843  0.15799077
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
#>            s1           s2         s3         s4          s5
#> gs1 0.1322743  0.122216769 0.06321831 0.15099473  0.18020272
#> gs2 0.1878729 -0.002001961 0.22666278 0.08462068 -0.08660604
#> gs3 0.1375218  0.131517876 0.13547157 0.11075140  0.08107227
#> gs4 0.1229140  0.169589360 0.10792946 0.14134695  0.23096335
#> gs5 0.1431006  0.009082125 0.09946726 0.21310843  0.15799077
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
#>      gsva 642.6850 645.5010 660.5564 660.8560 671.7195 689.8457    10
#>  bixverse 125.9492 126.3711 126.9739 126.7453 127.3844 128.3246    10
```
