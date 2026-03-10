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
#>                s1          s2          s3           s4           s5
#> gs1 -0.0638650428  0.07722360 -0.06203793 -0.006236476 -0.124861845
#> gs2  0.2212310494  0.18709054  0.10704192 -0.157679786  0.025438702
#> gs3  0.0108281822  0.07623144 -0.04423783 -0.164749297 -0.002397646
#> gs4  0.0009910496 -0.12307429  0.10897825 -0.068542190  0.132000696
#> gs5  0.0162243580 -0.24216373 -0.07598195 -0.057478720 -0.162065515
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>                s1          s2          s3          s4           s5
#> gs1 -0.0638753352  0.07722616 -0.06202476 -0.00623065 -0.124973921
#> gs2  0.2212310494  0.18707868  0.10702236 -0.15766665  0.025452841
#> gs3  0.0108352580  0.07622447 -0.04434669 -0.16475767 -0.002408061
#> gs4  0.0009910496 -0.12292891  0.10894471 -0.06843447  0.132003654
#> gs5  0.0162243580 -0.24216373 -0.07598195 -0.05747872 -0.162066617
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
#>      expr       min        lq      mean   median       uq       max neval
#>      gsva 2478.3308 2486.1406 2693.8697 2490.344 3144.306 3177.4011    10
#>  bixverse  570.0791  570.8198  572.6051  572.085  574.853  575.8698    10
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
#>             s1          s2            s3          s4            s5
#> gs1 -0.1147271 -0.08543100 -0.0386988776  0.07127380 -0.0486029435
#> gs2  0.3074678 -0.29125153  0.0821924633  0.02008761 -0.0316576809
#> gs3  0.3272081  0.21432306  0.0006195184 -0.07123876 -0.1734918157
#> gs4 -0.1916292 -0.08649112  0.1467308460  0.33845456 -0.1238031685
#> gs5 -0.6036840 -0.24508183  0.0493259381 -0.37402354  0.0002680626
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
#>             s1          s2            s3          s4            s5
#> gs1 -0.1147271 -0.08543100 -0.0386988776  0.07127380 -0.0486029435
#> gs2  0.3074678 -0.29125153  0.0821924633  0.02008761 -0.0316576809
#> gs3  0.3272081  0.21432306  0.0006195184 -0.07123876 -0.1734918157
#> gs4 -0.1916292 -0.08649112  0.1467308460  0.33845456 -0.1238031685
#> gs5 -0.6036840 -0.24508183  0.0493259381 -0.37402354  0.0002680626
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
#>      gsva 18.846641 18.863414 18.934912 18.870752 18.879030 19.552151    10
#>  bixverse  2.758365  2.760802  2.766917  2.764945  2.765814  2.791993    10
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
#>             s1          s2          s3         s4         s5
#> gs1 0.06413341  0.15118185  0.07009767 0.11997605 0.04629582
#> gs2 0.23169475  0.20116920  0.11974592 0.05641376 0.07204471
#> gs3 0.09400588  0.11688665  0.05253072 0.04983544 0.08959755
#> gs4 0.16555066  0.03992750  0.13784665 0.01444460 0.13531539
#> gs5 0.02852539 -0.01860348 -0.05944933 0.03674925 0.07897920
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
#>             s1          s2          s3         s4         s5
#> gs1 0.06413341  0.15118185  0.07009767 0.11997605 0.04629582
#> gs2 0.23169475  0.20116920  0.11974592 0.05641376 0.07204471
#> gs3 0.09400588  0.11688665  0.05253072 0.04983544 0.08959755
#> gs4 0.16555066  0.03992750  0.13784665 0.01444460 0.13531539
#> gs5 0.02852539 -0.01860348 -0.05944933 0.03674925 0.07897920
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
#>      gsva 614.8022 625.0392 702.1197 634.4265 656.8332 1308.6801    10
#>  bixverse 102.0209 102.5284 103.1024 103.0482 103.4619  104.2486    10
```
