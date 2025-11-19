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
#>              s1          s2            s3          s4          s5
#> gs1  0.11939563 -0.11811751 -0.2879038445 -0.03874867 -0.08273027
#> gs2  0.08291848 -0.37005508  0.4477334031  0.01543922  0.15755748
#> gs3 -0.28325673 -0.36485141 -0.1919007832  0.20333488  0.03406670
#> gs4  0.04761549 -0.07062547 -0.0002652648  0.06063840 -0.06983368
#> gs5 -0.14161068  0.21492048  0.1175303268 -0.12991197 -0.00809899
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1          s2            s3          s4           s5
#> gs1  0.11939423 -0.11824029 -0.2878953263 -0.03873681 -0.082726004
#> gs2  0.08291848 -0.37004751  0.4477334031  0.01544057  0.157676201
#> gs3 -0.28336294 -0.36485141 -0.1917990415  0.20344313  0.034064321
#> gs4  0.04762191 -0.07061325 -0.0002577195  0.06063711 -0.069833210
#> gs5 -0.14161682  0.21493924  0.1176538929 -0.12988840 -0.008081024
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
#>      expr      min        lq      mean    median        uq       max neval
#>      gsva 2478.060 2482.4354 2580.7739 2485.2592 2502.2976 3400.3131    10
#>  bixverse  499.016  499.3617  503.3965  500.5414  504.4689  520.8532    10
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
#>               s1          s2          s3          s4          s5
#> gs1 -0.008809846  0.07404341  0.08481443  0.09913239  0.01899850
#> gs2  0.060969112 -0.33724025 -0.15109017  0.04730359 -0.02529363
#> gs3  0.029652783 -0.33512994 -0.05034830 -0.08380039  0.26733032
#> gs4 -0.077992625 -0.02797281 -0.12474291 -0.08650759 -0.08575372
#> gs5  0.054265841  0.26697428  0.28321060  0.15091912 -0.31235534
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
#>               s1          s2          s3          s4          s5
#> gs1 -0.008809846  0.07404341  0.08481443  0.09913239  0.01899850
#> gs2  0.060969112 -0.33724025 -0.15109017  0.04730359 -0.02529363
#> gs3  0.029652783 -0.33512994 -0.05034830 -0.08380039  0.26733032
#> gs4 -0.077992625 -0.02797281 -0.12474291 -0.08650759 -0.08575372
#> gs5  0.054265841  0.26697428  0.28321060  0.15091912 -0.31235534
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
#>      expr       min        lq     mean    median        uq      max neval
#>      gsva 18.845234 18.854699 18.86284 18.858720 18.873312 18.88038    10
#>  bixverse  2.696595  2.700238  2.70279  2.702022  2.705687  2.70888    10
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
#>              s1          s2          s3         s4         s5
#> gs1  0.09529630  0.01226195 -0.03485719 0.11155915 0.08485097
#> gs2  0.16415781 -0.08978162  0.30392074 0.16152672 0.19438183
#> gs3  0.03026467 -0.04002123  0.04898387 0.20536208 0.10396828
#> gs4  0.12124818  0.01930436  0.04874923 0.19231468 0.04960162
#> gs5 -0.02899765  0.19030677  0.14698642 0.06077806 0.05766563
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
#>              s1          s2          s3         s4         s5
#> gs1  0.09529630  0.01226195 -0.03485719 0.11155915 0.08485097
#> gs2  0.16415781 -0.08978162  0.30392074 0.16152672 0.19438183
#> gs3  0.03026467 -0.04002123  0.04898387 0.20536208 0.10396828
#> gs4  0.12124818  0.01930436  0.04874923 0.19231468 0.04960162
#> gs5 -0.02899765  0.19030677  0.14698642 0.06077806 0.05766563
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
#>      gsva 635.0555 640.9978 649.8818 649.9616 653.7740 668.5987    10
#>  bixverse 117.6763 118.3727 118.6840 118.6133 119.1571 119.5048    10
```
