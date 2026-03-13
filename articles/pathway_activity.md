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
#> gs1  0.03686269  0.19494819 -0.11278927  0.18769563  0.05065022
#> gs2  0.03508263 -0.07545179 -0.26519825 -0.06310097  0.13823836
#> gs3 -0.18462912  0.03831025  0.03443377  0.10592651 -0.02201473
#> gs4  0.09461364 -0.07427210 -0.14171849 -0.06406094 -0.17101139
#> gs5 -0.19261159  0.15763791 -0.47230535 -0.15458954  0.10702596
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
#> gs1  0.03686762  0.19494394 -0.11278011  0.18768999  0.05064066
#> gs2  0.03509904 -0.07545167 -0.26510281 -0.06310736  0.13823300
#> gs3 -0.18463153  0.03831243  0.03445111  0.10591359 -0.02190520
#> gs4  0.09461364 -0.07427087 -0.14172571 -0.06405471 -0.17101358
#> gs5 -0.19260004  0.15761982 -0.47232361 -0.15456746  0.10692648
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
#>      gsva 2031.1409 2033.4653 2100.0177 2038.5905 2045.8135 2615.0265    10
#>  bixverse  640.9603  642.9439  646.6146  644.5229  646.2085  667.6438    10
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
#>               s1          s2          s3           s4          s5
#> gs1 -0.009359326  0.04345498  0.20730029 -0.023353317  0.01238404
#> gs2  0.128950657  0.15098043  0.14200464 -0.105503073 -0.17199108
#> gs3 -0.046277770 -0.05025991 -0.06227422 -0.090898408 -0.14131095
#> gs4  0.177275011  0.03300026  0.26403763  0.002285008 -0.16732278
#> gs5  0.326149210  0.07893128 -0.03024843  0.404225173  0.23783253
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
#>               s1          s2          s3           s4          s5
#> gs1 -0.009359326  0.04345498  0.20730029 -0.023353317  0.01238404
#> gs2  0.128950657  0.15098043  0.14200464 -0.105503073 -0.17199108
#> gs3 -0.046277770 -0.05025991 -0.06227422 -0.090898408 -0.14131095
#> gs4  0.177275011  0.03300026  0.26403763  0.002285008 -0.16732278
#> gs5  0.326149210  0.07893128 -0.03024843  0.404225173  0.23783253
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
#>      expr       min        lq      mean   median       uq      max neval
#>      gsva 16.381409 16.388887 16.527064 16.41507 16.42617 17.01935    10
#>  bixverse  2.839237  2.841619  2.843473  2.84314  2.84513  2.84782    10
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
#>             s1         s2          s3         s4           s5
#> gs1 0.10012020 0.21648786  0.02513695 0.20472021  0.102343464
#> gs2 0.16262722 0.08272948 -0.04443606 0.07332105  0.147577189
#> gs3 0.06535876 0.13842006  0.13049957 0.12576408  0.079610277
#> gs4 0.15468338 0.10867078  0.04692916 0.04987452 -0.002227429
#> gs5 0.02905189 0.15690685 -0.20547742 0.07104171  0.155643858
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
#>             s1         s2          s3         s4           s5
#> gs1 0.10012020 0.21648786  0.02513695 0.20472021  0.102343464
#> gs2 0.16262722 0.08272948 -0.04443606 0.07332105  0.147577189
#> gs3 0.06535876 0.13842006  0.13049957 0.12576408  0.079610277
#> gs4 0.15468338 0.10867078  0.04692916 0.04987452 -0.002227429
#> gs5 0.02905189 0.15690685 -0.20547742 0.07104171  0.155643858
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
#>      expr      min       lq     mean   median       uq       max neval
#>      gsva 604.5376 611.9876 739.9460 625.4833 635.2711 1224.5165    10
#>  bixverse 129.7344 130.4624 133.5867 131.3484 135.0119  147.6294    10
```
