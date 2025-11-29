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
#>              s1           s2          s3          s4           s5
#> gs1  0.42105572 -0.060707309 -0.53338740 -0.49198353 -0.356661546
#> gs2 -0.07445493 -0.017467731 -0.08787027  0.18305864  0.004608677
#> gs3 -0.15239708 -0.007203049  0.23983295  0.43043903 -0.058201981
#> gs4  0.31255211 -0.235247872 -0.01267029  0.04954988 -0.164434181
#> gs5 -0.07130246 -0.066562441 -0.02169517  0.12541076  0.101594334
```

If we compare against the original implementation of GSVA

``` r
library(GSVA)

gsvaPar <- gsvaParam(X, gs)

gsva_res_gaussian <- as.matrix(
  gsva(gsvaPar, verbose=FALSE)
)

gsva_res_gaussian[1:5, 1:5]
#>              s1           s2          s3          s4           s5
#> gs1  0.42108195 -0.060741482 -0.53338740 -0.49201765 -0.356640663
#> gs2 -0.07446244 -0.017577942 -0.08784127  0.18305677  0.004607153
#> gs3 -0.15251468 -0.007199846  0.23983295  0.43044919 -0.058197828
#> gs4  0.31252704 -0.235273714 -0.01265327  0.04955974 -0.164434181
#> gs5 -0.07130275 -0.066555265 -0.02179375  0.12541342  0.101594665
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
#>      gsva 2457.1831 2464.6465 2553.3194 2465.7472 2479.1585 3305.8034    10
#>  bixverse  480.4329  480.7238  481.9848  481.6229  482.8805  485.6477    10
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
#> gs1  0.36187136 -0.01496766 -0.08064324  0.12477001  0.39107556
#> gs2  0.01689041  0.13168303  0.19614987  0.31048310  0.18382534
#> gs3 -0.10338828 -0.28823323 -0.02908575  0.15263298 -0.02621127
#> gs4 -0.18726507 -0.05574603 -0.32223292 -0.01938698 -0.05570599
#> gs5  0.01071131  0.19335855 -0.03194207 -0.18629991 -0.20819117
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
#> gs1  0.36187136 -0.01496766 -0.08064324  0.12477001  0.39107556
#> gs2  0.01689041  0.13168303  0.19614987  0.31048310  0.18382534
#> gs3 -0.10338828 -0.28823323 -0.02908575  0.15263298 -0.02621127
#> gs4 -0.18726507 -0.05574603 -0.32223292 -0.01938698 -0.05570599
#> gs5  0.01071131  0.19335855 -0.03194207 -0.18629991 -0.20819117
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
#>      expr       min        lq      mean    median        uq      max neval
#>      gsva 18.831555 18.835555 18.842753 18.837974 18.850031 18.86321    10
#>  bixverse  2.702679  2.704643  2.710554  2.706018  2.706531  2.73660    10
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
#> gs1 0.35669206  0.17397610 -0.08946712 -0.1567727 -0.10362544
#> gs2 0.08394799  0.11231188  0.10126876  0.2113155  0.10898206
#> gs3 0.11228037  0.19191815  0.20685201  0.3052967  0.03767693
#> gs4 0.24515689 -0.06885469  0.11998731  0.1237157 -0.03535827
#> gs5 0.09299553  0.07519045  0.11526918  0.1729491  0.12896766
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
#> gs1 0.35669206  0.17397610 -0.08946712 -0.1567727 -0.10362544
#> gs2 0.08394799  0.11231188  0.10126876  0.2113155  0.10898206
#> gs3 0.11228037  0.19191815  0.20685201  0.3052967  0.03767693
#> gs4 0.24515689 -0.06885469  0.11998731  0.1237157 -0.03535827
#> gs5 0.09299553  0.07519045  0.11526918  0.1729491  0.12896766
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
#>      gsva 603.9823 615.6363 620.6967 621.3493 627.5297 636.0291    10
#>  bixverse 111.0514 111.1607 112.1578 112.4541 112.9089 113.2650    10
```
