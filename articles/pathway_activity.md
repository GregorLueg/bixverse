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
#> gs1 -0.17044896 -0.06266398  0.05343558  0.06700869 -0.03515929
#> gs2  0.14360906 -0.13156188 -0.03288786  0.09583991 -0.14195712
#> gs3 -0.10671876 -0.08644815  0.00345292  0.02007892  0.08499326
#> gs4 -0.02658630  0.05744306 -0.19427813 -0.10539521  0.25465283
#> gs5 -0.07879267  0.09742901 -0.04279042  0.17162354 -0.12033657
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
#> gs1 -0.17045413 -0.06265001  0.053330099  0.06700572 -0.03505343
#> gs2  0.14360906 -0.13143331 -0.032894209  0.09584871 -0.14195635
#> gs3 -0.10672636 -0.08645214  0.003453436  0.02007010  0.08499621
#> gs4 -0.02658425  0.05746244 -0.194264731 -0.10540757  0.25464311
#> gs5 -0.07880687  0.09742339 -0.042782910  0.17163637 -0.12034533
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
#>      gsva 2467.891 2476.6241 2569.3375 2479.9338 2488.1521 3348.9284    10
#>  bixverse  498.093  498.7967  500.1011  499.2617  499.6569  505.8774    10
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
#> gs1 -0.19549004  0.07537397 -0.18269444  0.05716237  0.13480901
#> gs2  0.02239211  0.24130960  0.19261994 -0.15235948 -0.02496863
#> gs3  0.01480581 -0.09434441  0.08571427  0.01779511 -0.11429316
#> gs4  0.19520412 -0.17576493 -0.14649328 -0.38527189  0.19101752
#> gs5  0.09631913  0.04426459 -0.06609087  0.11195183 -0.18101182
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
#> gs1 -0.19549004  0.07537397 -0.18269444  0.05716237  0.13480901
#> gs2  0.02239211  0.24130960  0.19261994 -0.15235948 -0.02496863
#> gs3  0.01480581 -0.09434441  0.08571427  0.01779511 -0.11429316
#> gs4  0.19520412 -0.17576493 -0.14649328 -0.38527189  0.19101752
#> gs5  0.09631913  0.04426459 -0.06609087  0.11195183 -0.18101182
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
#>      gsva 18.850206 18.857242 18.864206 18.860412 18.868831 18.888242    10
#>  bixverse  2.713348  2.719535  2.719907  2.720284  2.722573  2.724803    10
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
#>                s1         s2          s3         s4         s5
#> gs1 -0.0242072032 0.07576582  0.12708650 0.16669725 0.09589627
#> gs2  0.1789042545 0.03075465  0.05238697 0.11558092 0.03667857
#> gs3  0.0478415331 0.09371521  0.13515926 0.08072932 0.14183171
#> gs4  0.0791247558 0.14352720 -0.01487278 0.07755102 0.30854538
#> gs5  0.0005355318 0.11387354  0.02925914 0.23145338 0.06985116
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
#>                s1         s2          s3         s4         s5
#> gs1 -0.0242072032 0.07576582  0.12708650 0.16669725 0.09589627
#> gs2  0.1789042545 0.03075465  0.05238697 0.11558092 0.03667857
#> gs3  0.0478415331 0.09371521  0.13515926 0.08072932 0.14183171
#> gs4  0.0791247558 0.14352720 -0.01487278 0.07755102 0.30854538
#> gs5  0.0005355318 0.11387354  0.02925914 0.23145338 0.06985116
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
#>      expr      min       lq     mean  median       uq      max neval
#>      gsva 658.3453 680.9555 683.5526 683.269 693.3116 697.7402    10
#>  bixverse 103.2901 103.9918 104.5418 104.645 105.1309 105.6338    10
```
