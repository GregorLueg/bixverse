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
#>              s1           s2            s3           s4           s5
#> gs1 -0.09480484 -0.270244298  4.379277e-05 -0.226720256  0.051010674
#> gs2 -0.06872688 -0.185850304 -1.752628e-01  0.163209719 -0.044254153
#> gs3 -0.17039256  0.002742701 -4.798406e-01  0.008612669 -0.143951134
#> gs4 -0.05967178 -0.022660054 -5.376422e-02 -0.084452701  0.002624792
#> gs5 -0.16296488  0.093331681  2.402743e-01 -0.049245219 -0.100997505
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1           s2            s3           s4           s5
#> gs1 -0.09480055 -0.270227723  3.346876e-05 -0.226718676  0.051014363
#> gs2 -0.06872820 -0.185859451 -1.752512e-01  0.163221127 -0.044253296
#> gs3 -0.17037212  0.002742701 -4.798591e-01  0.008612669 -0.144126340
#> gs4 -0.05967927 -0.022654984 -5.366943e-02 -0.084453558  0.002621099
#> gs5 -0.16296488  0.093343530  2.402743e-01 -0.049223372 -0.100992711
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
#>      gsva 2503.9109 2509.2728 2524.9728 2517.6866 2534.8203 2581.6558    10
#>  bixverse  481.8508  482.5793  483.1429  483.2812  483.4722  484.3998    10
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
#> gs1 -0.08976593 -0.03411952  0.06763712 -0.15913272  0.01575779
#> gs2 -0.12221290  0.07336711 -0.06816841  0.16229418 -0.17011497
#> gs3 -0.10679174  0.16138863 -0.28242331  0.36922636 -0.06118072
#> gs4  0.28293548  0.09674616 -0.00846292 -0.09243402  0.08382052
#> gs5 -0.18220566  0.03738355  0.02203415 -0.10652119 -0.08567389
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
#> gs1 -0.08976593 -0.03411952  0.06763712 -0.15913272  0.01575779
#> gs2 -0.12221290  0.07336711 -0.06816841  0.16229418 -0.17011497
#> gs3 -0.10679174  0.16138863 -0.28242331  0.36922636 -0.06118072
#> gs4  0.28293548  0.09674616 -0.00846292 -0.09243402  0.08382052
#> gs5 -0.18220566  0.03738355  0.02203415 -0.10652119 -0.08567389
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
#>      gsva 18.868157 18.871875 18.96718 18.889263 18.892764 19.70879    10
#>  bixverse  2.684808  2.688096  2.69192  2.689882  2.693639  2.71309    10
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
#>               s1          s2          s3           s4          s5
#> gs1  0.055923481 -0.03546097  0.14274910 -0.002263917  0.16909662
#> gs2  0.070602124  0.04822221  0.04685610  0.184340119  0.10341196
#> gs3 -0.039171780  0.05472790 -0.22142726  0.119385456 -0.03452588
#> gs4  0.012802014  0.08520621  0.03806565  0.134075055  0.07110518
#> gs5  0.004406202  0.14612637  0.20857382 -0.004540797  0.04230400
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
#>               s1          s2          s3           s4          s5
#> gs1  0.055923481 -0.03546097  0.14274910 -0.002263917  0.16909662
#> gs2  0.070602124  0.04822221  0.04685610  0.184340119  0.10341196
#> gs3 -0.039171780  0.05472790 -0.22142726  0.119385456 -0.03452588
#> gs4  0.012802014  0.08520621  0.03806565  0.134075055  0.07110518
#> gs5  0.004406202  0.14612637  0.20857382 -0.004540797  0.04230400
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
#>      gsva 609.4235 626.9762 811.7209 638.1435 651.5413 1534.7625    10
#>  bixverse 102.3169 102.5287 103.2069 102.9441 103.1264  106.6607    10
```
