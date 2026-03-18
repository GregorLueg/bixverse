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
#> gs1  0.07680500 -0.11302533 -0.03934923  0.08032242  0.22772860
#> gs2  0.06230482  0.01609217  0.04576508 -0.12423160  0.21910595
#> gs3 -0.13359384  0.16966518 -0.07890620 -0.04965710 -0.03368295
#> gs4  0.33888260 -0.12072944 -0.09695866  0.05321099 -0.16570657
#> gs5 -0.02333842  0.01419487  0.13741696  0.14782561 -0.09745488
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
#> gs1  0.07679513 -0.11304124 -0.03935872  0.08031931  0.22771828
#> gs2  0.06231932  0.01609697  0.04575791 -0.12421915  0.21910286
#> gs3 -0.13359060  0.16967018 -0.07891350 -0.04966425 -0.03367729
#> gs4  0.33881438 -0.12071677 -0.09704694  0.05316119 -0.16574295
#> gs5 -0.02341718  0.01417917  0.13741696  0.14781638 -0.09745179
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
  times = 5L
)
#> Unit: milliseconds
#>      expr       min        lq      mean   median        uq       max neval
#>      gsva 2485.5005 2497.8297 2506.0742 2501.559 2503.1027 2542.3790     5
#>  bixverse  570.5823  573.0117  573.3832  573.490  574.1456  575.6863     5
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
#>              s1          s2          s3          s4            s5
#> gs1  0.03555496 -0.09479711 -0.10109203 -0.02258372 -0.0009649248
#> gs2 -0.05455040  0.08703022 -0.03047756 -0.03334293  0.0902819622
#> gs3 -0.09073528 -0.19721119 -0.06507774 -0.07967113 -0.0565364066
#> gs4 -0.11190642  0.25317653  0.01823319 -0.08895262  0.2206121197
#> gs5  0.18177513 -0.11916659  0.25600796  0.32930743 -0.3144384229
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
#>              s1          s2          s3          s4            s5
#> gs1  0.03555496 -0.09479711 -0.10109203 -0.02258372 -0.0009649248
#> gs2 -0.05455040  0.08703022 -0.03047756 -0.03334293  0.0902819622
#> gs3 -0.09073528 -0.19721119 -0.06507774 -0.07967113 -0.0565364066
#> gs4 -0.11190642  0.25317653  0.01823319 -0.08895262  0.2206121197
#> gs5  0.18177513 -0.11916659  0.25600796  0.32930743 -0.3144384229
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
  times = 5L
)
#> Unit: seconds
#>      expr      min        lq      mean    median        uq       max neval
#>      gsva 18.85632 18.880832 18.885536 18.890584 18.892013 18.907932     5
#>  bixverse  2.77629  2.777754  2.778904  2.778397  2.778452  2.783629     5
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
#>             s1          s2          s3         s4          s5
#> gs1 0.16757295  0.02291227  0.05352661 0.15503778  0.17662940
#> gs2 0.15706902  0.12685854  0.17477382 0.04851582  0.15316529
#> gs3 0.04328484  0.18521995  0.03643323 0.07810993  0.12536501
#> gs4 0.20389567 -0.11191251 -0.04096866 0.16114398  0.06374350
#> gs5 0.13401626  0.14820729  0.12236692 0.13516653 -0.01950481
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
#>             s1          s2          s3         s4          s5
#> gs1 0.16757295  0.02291227  0.05352661 0.15503778  0.17662940
#> gs2 0.15706902  0.12685854  0.17477382 0.04851582  0.15316529
#> gs3 0.04328484  0.18521995  0.03643323 0.07810993  0.12536501
#> gs4 0.20389567 -0.11191251 -0.04096866 0.16114398  0.06374350
#> gs5 0.13401626  0.14820729  0.12236692 0.13516653 -0.01950481
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
  times = 5L
)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>      gsva 615.9093 619.9724 625.3839 624.7308 624.9001 641.4067     5
#>  bixverse 121.9735 122.4048 126.7067 122.7186 123.1249 143.3114     5
```
