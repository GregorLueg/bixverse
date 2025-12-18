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
#>              s1          s2          s3          s4          s5
#> gs1 -0.16217577 -0.07518150 -0.03987554  0.01792652 -0.04335536
#> gs2 -0.25051360  0.27700731  0.04050408  0.10224711 -0.23343932
#> gs3  0.10206806  0.08516579 -0.07634456 -0.20059948 -0.20669210
#> gs4 -0.09333005  0.23367194 -0.07284775 -0.26152567 -0.17032866
#> gs5  0.04673594  0.18112160 -0.39819724 -0.06408744  0.01209156
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1          s2          s3          s4          s5
#> gs1 -0.16207777 -0.07506555 -0.03972891  0.01792701 -0.04346236
#> gs2 -0.25051951  0.27701838  0.04049968  0.10224711 -0.23344805
#> gs3  0.10205214  0.08505883 -0.07633478 -0.20070195 -0.20680205
#> gs4 -0.09343794  0.23367626 -0.07285915 -0.26163409 -0.17032597
#> gs5  0.04671056  0.18112160 -0.39819724 -0.06409101  0.01196597
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
#>      expr      min       lq     mean    median        uq       max neval
#>      gsva 2488.906 2495.739 2503.223 2498.8739 2502.6135 2544.8775    10
#>  bixverse  481.068  481.466  488.694  482.6694  483.7827  544.7485    10
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
#>              s1          s2           s3           s4         s5
#> gs1 -0.03932111  0.01494481 -0.134320891 -0.008100723 -0.3047451
#> gs2  0.06704029 -0.03497663  0.222670285 -0.070530404 -0.2427219
#> gs3 -0.04878407  0.04151562 -0.005210885 -0.053318387 -0.0517624
#> gs4  0.33908178 -0.03068597  0.160271574  0.280315869  0.0466641
#> gs5 -0.11122965 -0.22048840  0.034984260 -0.050050844 -0.2343698
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
#>              s1          s2           s3           s4         s5
#> gs1 -0.03932111  0.01494481 -0.134320891 -0.008100723 -0.3047451
#> gs2  0.06704029 -0.03497663  0.222670285 -0.070530404 -0.2427219
#> gs3 -0.04878407  0.04151562 -0.005210885 -0.053318387 -0.0517624
#> gs4  0.33908178 -0.03068597  0.160271574  0.280315869  0.0466641
#> gs5 -0.11122965 -0.22048840  0.034984260 -0.050050844 -0.2343698
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
#>      expr       min        lq      mean    median        uq       max neval
#>      gsva 18.840568 18.845550 18.935754 18.854842 18.860833 19.684526    10
#>  bixverse  2.677031  2.678379  2.681148  2.680817  2.684341  2.684947    10
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
#>             s1        s2          s3          s4         s5
#> gs1 0.03807339 0.0419132  0.11258750  0.10559899 0.07411777
#> gs2 0.02783148 0.2110984  0.14247690  0.19332503 0.04974275
#> gs3 0.14504884 0.1316739  0.06717493  0.03719170 0.01848772
#> gs4 0.04390467 0.1662134  0.05322878 -0.03532138 0.02651937
#> gs5 0.03171569 0.1560405 -0.09680967  0.09950538 0.14837551
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
#>             s1        s2          s3          s4         s5
#> gs1 0.03807339 0.0419132  0.11258750  0.10559899 0.07411777
#> gs2 0.02783148 0.2110984  0.14247690  0.19332503 0.04974275
#> gs3 0.14504884 0.1316739  0.06717493  0.03719170 0.01848772
#> gs4 0.04390467 0.1662134  0.05322878 -0.03532138 0.02651937
#> gs5 0.03171569 0.1560405 -0.09680967  0.09950538 0.14837551
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
#>      expr       min        lq     mean   median       uq       max neval
#>      gsva 601.97345 619.85635 719.9089 640.0571 672.5145 1432.7473    10
#>  bixverse  97.90607  99.36771 100.2404 100.2935 100.6506  103.3355    10
```
