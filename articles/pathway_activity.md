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
n <- 100 # number of samples

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
#>              s1           s2          s3          s4         s5
#> gs1 -0.01320378 -0.093492292  0.10111879 -0.10211421  0.1019284
#> gs2  0.14308516 -0.132894778  0.13513867  0.09514166 -0.1564986
#> gs3 -0.21491357 -0.075694242  0.12316338 -0.04042700 -0.1765931
#> gs4  0.12565312  0.129707548 -0.07585571 -0.22496932  0.2923787
#> gs5  0.15323107  0.002392075 -0.07922392 -0.09677161  0.1677735
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose = FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1           s2          s3          s4         s5
#> gs1 -0.01320378 -0.093486863  0.10113341 -0.10211972  0.1019327
#> gs2  0.14309364 -0.133006505  0.13515390  0.09503537 -0.1563926
#> gs3 -0.21491357 -0.075675821  0.12306305 -0.04043898 -0.1766003
#> gs4  0.12565056  0.129698884 -0.07586936 -0.22497445  0.2923787
#> gs5  0.15323059  0.002390686 -0.07921143 -0.09676413  0.1677762
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
  times = 3L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq      max neval
#>      gsva 2505.7374 2511.5945 2525.2353 2517.4516 2534.9843 2552.517     3
#>  bixverse  580.6615  580.7829  581.0223  580.9043  581.2026  581.501     3
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
#>              s1          s2           s3          s4           s5
#> gs1 -0.13315522  0.11296689  0.028935180 -0.03946902  0.009183941
#> gs2  0.02985661 -0.04947534  0.133505116  0.05730318  0.089421030
#> gs3  0.06049196 -0.09175719 -0.097897108  0.15923866 -0.227262363
#> gs4  0.12682541  0.00356020 -0.024443641 -0.04817659  0.332142450
#> gs5  0.04596053 -0.13591350 -0.008891869  0.01180001 -0.117298668
```

This is how you would do this in the official GSVA code:

``` r
# update this to use for kcdf = "Poisson"
gsvaParPoisson <- gsvaParam(X1, gs, kcdf = "Poisson")

gsva_res_poisson <- as.matrix(
  gsva(gsvaParPoisson, verbose = FALSE)
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
#>              s1          s2           s3          s4           s5
#> gs1 -0.13315522  0.11296689  0.028935180 -0.03946902  0.009183941
#> gs2  0.02985661 -0.04947534  0.133505116  0.05730318  0.089421030
#> gs3  0.06049196 -0.09175719 -0.097897108  0.15923866 -0.227262363
#> gs4  0.12682541  0.00356020 -0.024443641 -0.04817659  0.332142450
#> gs5  0.04596053 -0.13591350 -0.008891869  0.01180001 -0.117298668
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
  times = 3L
)
#> Unit: seconds
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 18.893905 18.896677 18.899481 18.899448 18.902270 18.905091     3
#>  bixverse  2.766826  2.771615  2.777271  2.776405  2.782494  2.788582     3
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
#>              s1         s2         s3          s4         s5
#> gs1  0.11199829 0.06760971 0.15310242  0.06196397 0.15278275
#> gs2  0.12984610 0.07205401 0.13272593  0.11618010 0.04902417
#> gs3 -0.03134935 0.03046895 0.21958267  0.03765917 0.06731794
#> gs4  0.25277204 0.14904210 0.09270736 -0.03693458 0.31172069
#> gs5  0.18062070 0.09837657 0.10226841  0.06596902 0.19550864
```

This is the way you would run the ssGSEA algorithm in the GSVA package.

``` r
# update this to the ssgsea parameters
ssgseaPar <- ssgseaParam(X, gs)

ssgsea_res <- as.matrix(
  gsva(ssgseaPar, verbose = FALSE)
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
#>              s1         s2         s3          s4         s5
#> gs1  0.11199829 0.06760971 0.15310242  0.06196397 0.15278275
#> gs2  0.12984610 0.07205401 0.13272593  0.11618010 0.04902417
#> gs3 -0.03134935 0.03046895 0.21958267  0.03765917 0.06731794
#> gs4  0.25277204 0.14904210 0.09270736 -0.03693458 0.31172069
#> gs5  0.18062070 0.09837657 0.10226841  0.06596902 0.19550864
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
  times = 3L
)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>      gsva 620.2673 623.9776 625.2376 627.6879 627.7227 627.7575     3
#>  bixverse 119.0858 119.9084 120.4692 120.7310 121.1609 121.5908     3
```
