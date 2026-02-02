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
#>               s1           s2          s3           s4           s5
#> gs1 -0.179783193 -0.324056403 0.002646492  0.135305352 -0.007232576
#> gs2 -0.129740090  0.001673989 0.013618344 -0.126805769  0.052353083
#> gs3  0.273429192  0.004027738 0.102973077  0.427661840  0.178132704
#> gs4 -0.009012528  0.002005581 0.083458489  0.002266909 -0.011771512
#> gs5 -0.007770160  0.283630073 0.104347408 -0.266056617  0.111906039
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>               s1           s2         s3           s4           s5
#> gs1 -0.179888856 -0.324061084 0.00266080  0.135310388 -0.007225557
#> gs2 -0.129745359  0.001769447 0.01351016 -0.126793131  0.052354049
#> gs3  0.273419506  0.003879021 0.10309282  0.427780617  0.178132704
#> gs4 -0.009012455  0.002001776 0.08346251  0.002270824 -0.011762540
#> gs5 -0.007794797  0.283739772 0.10436321 -0.266071966  0.111885506
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
#>      expr       min        lq     mean    median       uq       max neval
#>      gsva 2525.7687 2537.7017 2841.446 2545.9676 3529.105 3544.5797    10
#>  bixverse  572.2308  573.0029  574.229  574.2771  574.852  576.9246    10
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
#>               s1          s2          s3          s4           s5
#> gs1  0.006486357 -0.25624310 -0.04865821  0.04668242  0.206267341
#> gs2 -0.113884127  0.16313245  0.01566856  0.22329577  0.215286461
#> gs3  0.189586487 -0.03610791  0.15381537 -0.38895489 -0.189826228
#> gs4  0.112920026 -0.12546301  0.16778171 -0.16108323 -0.122970622
#> gs5  0.018052942 -0.21385294  0.06717491 -0.18588812  0.007938484
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
#>               s1          s2          s3          s4           s5
#> gs1  0.006486357 -0.25624310 -0.04865821  0.04668242  0.206267341
#> gs2 -0.113884127  0.16313245  0.01566856  0.22329577  0.215286461
#> gs3  0.189586487 -0.03610791  0.15381537 -0.38895489 -0.189826228
#> gs4  0.112920026 -0.12546301  0.16778171 -0.16108323 -0.122970622
#> gs5  0.018052942 -0.21385294  0.06717491 -0.18588812  0.007938484
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
#>      gsva 18.897255 18.903852 19.100212 18.921462 18.996759 19.801889    10
#>  bixverse  2.766184  2.770214  2.772946  2.771573  2.773226  2.790226    10
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
#>             s1          s2         s3          s4         s5
#> gs1 0.04660481 -0.08611041 0.09827951  0.14990469 0.11856022
#> gs2 0.02668100  0.08695356 0.08914728  0.03283895 0.07026851
#> gs3 0.30621798  0.24817764 0.17654645  0.29231028 0.16952853
#> gs4 0.11242266  0.14253486 0.14029347  0.12446179 0.11033051
#> gs5 0.17168248  0.22867215 0.11828895 -0.03563434 0.17890117
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
#>             s1          s2         s3          s4         s5
#> gs1 0.04660481 -0.08611041 0.09827951  0.14990469 0.11856022
#> gs2 0.02668100  0.08695356 0.08914728  0.03283895 0.07026851
#> gs3 0.30621798  0.24817764 0.17654645  0.29231028 0.16952853
#> gs4 0.11242266  0.14253486 0.14029347  0.12446179 0.11033051
#> gs5 0.17168248  0.22867215 0.11828895 -0.03563434 0.17890117
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
#>      gsva 638.1720 649.4125 660.2691 656.8052 665.9964 693.2247    10
#>  bixverse 102.2992 103.0438 103.7283 103.4976 103.7496 106.2961    10
```
