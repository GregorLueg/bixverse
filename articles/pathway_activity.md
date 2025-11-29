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
#>               s1          s2          s3         s4           s5
#> gs1 -0.059587213  0.15066269  0.15785470 0.28629161 -0.002956328
#> gs2  0.130636123 -0.02780410 -0.01144541 0.13454750 -0.236769046
#> gs3  0.012890764 -0.10082289  0.15017491 0.24751645 -0.152289071
#> gs4 -0.052736167  0.03885136  0.09501237 0.11677047 -0.091179977
#> gs5  0.001186552 -0.24184258 -0.04222997 0.06673762  0.223143282
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>               s1          s2          s3         s4          s5
#> gs1 -0.059595089  0.15065010  0.15785655 0.28628677 -0.00295354
#> gs2  0.130634065 -0.02778239 -0.01137134 0.13456425 -0.23673395
#> gs3  0.012890764 -0.10086329  0.15015453 0.24754307 -0.15252540
#> gs4 -0.052746241  0.03883864  0.09502088 0.11676030 -0.09118359
#> gs5  0.001180693 -0.24184362 -0.04222978 0.06675151  0.22304152
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
#>      expr       min        lq      mean   median        uq       max neval
#>      gsva 2465.2730 2470.4666 2561.2095 2472.795 2474.6403 3332.5254    10
#>  bixverse  472.0914  472.5912  473.0612  472.898  473.7762  474.3378    10
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
#>              s1          s2           s3          s4          s5
#> gs1  0.12015954  0.12661374 -0.002463926 -0.25065544 -0.07621185
#> gs2 -0.37542403  0.01973761 -0.426346384  0.11443177  0.26001611
#> gs3 -0.09414836 -0.44837370  0.056631610  0.22071040  0.09082449
#> gs4 -0.21320436  0.03404501 -0.067371829 -0.02843437  0.19630530
#> gs5 -0.09081271 -0.13053859  0.157610998 -0.04123752 -0.10097780
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
#>              s1          s2           s3          s4          s5
#> gs1  0.12015954  0.12661374 -0.002463926 -0.25065544 -0.07621185
#> gs2 -0.37542403  0.01973761 -0.426346384  0.11443177  0.26001611
#> gs3 -0.09414836 -0.44837370  0.056631610  0.22071040  0.09082449
#> gs4 -0.21320436  0.03404501 -0.067371829 -0.02843437  0.19630530
#> gs5 -0.09081271 -0.13053859  0.157610998 -0.04123752 -0.10097780
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
#>      expr      min        lq      mean   median        uq       max neval
#>      gsva 18.84520 18.851476 18.859687 18.85897 18.868816 18.874153    10
#>  bixverse  2.68475  2.686069  2.687171  2.68680  2.688289  2.691155    10
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
#>             s1           s2          s3        s4          s5
#> gs1 0.05926982  0.187105610  0.25105325 0.2390465  0.09583468
#> gs2 0.15467666  0.133515592 -0.01216673 0.1121622 -0.11206098
#> gs3 0.11489652  0.115054692  0.13845078 0.2869225  0.07663045
#> gs4 0.07459662  0.165794593  0.18585565 0.1565263  0.05346765
#> gs5 0.13469942 -0.009127753  0.10222149 0.1709920  0.18153912
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
#>             s1           s2          s3        s4          s5
#> gs1 0.05926982  0.187105610  0.25105325 0.2390465  0.09583468
#> gs2 0.15467666  0.133515592 -0.01216673 0.1121622 -0.11206098
#> gs3 0.11489652  0.115054692  0.13845078 0.2869225  0.07663045
#> gs4 0.07459662  0.165794593  0.18585565 0.1565263  0.05346765
#> gs5 0.13469942 -0.009127753  0.10222149 0.1709920  0.18153912
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
#>      gsva 625.1906 626.3972 635.7114 630.0289 645.8510 658.2132    10
#>  bixverse 129.8343 129.8866 131.0622 130.7722 131.9632 133.5772    10
```
