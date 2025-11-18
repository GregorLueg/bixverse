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
#>              s1          s2           s3          s4          s5
#> gs1 -0.10414096  0.06198255 -0.131958977 -0.11219463  0.10426769
#> gs2 -0.06256114 -0.15974090 -0.013662502  0.06314011  0.08135172
#> gs3 -0.04555765 -0.28084360 -0.122844729  0.02967194  0.30373554
#> gs4 -0.25708916 -0.11252141  0.069027782 -0.28074188  0.03496737
#> gs5 -0.23058215 -0.19703249 -0.005814211  0.21048751 -0.33644934
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1          s2           s3          s4          s5
#> gs1 -0.10402687  0.06197561 -0.131957479 -0.11240382  0.10428047
#> gs2 -0.06256771 -0.15973555 -0.013675656  0.06314449  0.08136305
#> gs3 -0.04564959 -0.28092993 -0.122844729  0.02967194  0.30374133
#> gs4 -0.25708446 -0.11251812  0.069033053 -0.28061828  0.03497937
#> gs5 -0.23039629 -0.19712626 -0.005796247  0.21059716 -0.33644934
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
#>      gsva 2462.2631 2465.5159 2561.5061 2470.4733 2490.4844 3325.1976    10
#>  bixverse  497.3681  497.9669  498.9191  499.0746  500.0164  500.1218    10
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
#> gs1 -0.03732752 -0.19871575 -0.02172686  0.14817513  0.04360175
#> gs2  0.12082599 -0.07621550  0.20116777  0.08858233  0.12954411
#> gs3  0.17306534  0.21404369 -0.13691868 -0.09559570 -0.10596296
#> gs4 -0.13270192  0.06885957  0.03848889  0.22493063 -0.20856231
#> gs5 -0.24364003 -0.21362195 -0.34754850  0.22844097 -0.02121493
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
#> gs1 -0.03732752 -0.19871575 -0.02172686  0.14817513  0.04360175
#> gs2  0.12082599 -0.07621550  0.20116777  0.08858233  0.12954411
#> gs3  0.17306534  0.21404369 -0.13691868 -0.09559570 -0.10596296
#> gs4 -0.13270192  0.06885957  0.03848889  0.22493063 -0.20856231
#> gs5 -0.24364003 -0.21362195 -0.34754850  0.22844097 -0.02121493
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
#>      gsva 18.82741 18.831300 18.843885 18.834288 18.853115 18.877464    10
#>  bixverse  2.69714  2.698494  2.701179  2.701832  2.703123  2.705856    10
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
#>               s1          s2         s3          s4          s5
#> gs1 0.0306837337  0.20934833 0.05937760  0.06995628 0.165018323
#> gs2 0.0410987761  0.09596220 0.14040351  0.20988164 0.113395840
#> gs3 0.1838608428 -0.05956005 0.21127686  0.17431571 0.378124617
#> gs4 0.0311369032  0.06684633 0.14462451 -0.02621928 0.099702947
#> gs5 0.0007043416  0.03350229 0.08687154  0.17149717 0.004852471
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
#>               s1          s2         s3          s4          s5
#> gs1 0.0306837337  0.20934833 0.05937760  0.06995628 0.165018323
#> gs2 0.0410987761  0.09596220 0.14040351  0.20988164 0.113395840
#> gs3 0.1838608428 -0.05956005 0.21127686  0.17431571 0.378124617
#> gs4 0.0311369032  0.06684633 0.14462451 -0.02621928 0.099702947
#> gs5 0.0007043416  0.03350229 0.08687154  0.17149717 0.004852471
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
#>      expr       min       lq     mean   median       uq      max neval
#>      gsva 615.26983 619.3168 626.4831 624.2660 631.5253 648.1907    10
#>  bixverse  99.68772 100.0765 100.8924 100.5515 101.1318 103.9173    10
```
