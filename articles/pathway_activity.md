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
#> gs1  0.07138880  0.06098404 -0.08629830 -0.17949432  0.07135021
#> gs2  0.05245457 -0.11197691 -0.14411013 -0.20334903  0.01451552
#> gs3 -0.19453865  0.21379668 -0.18038804 -0.13957182 -0.19982854
#> gs4  0.31203690 -0.02125201 -0.25247325  0.02559055 -0.33436258
#> gs5  0.18774145 -0.06790280 -0.02898929 -0.08483732  0.09957309
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
#> gs1  0.07138667  0.06099006 -0.08630301 -0.17949886  0.07145038
#> gs2  0.05246474 -0.11197691 -0.14411679 -0.20334477  0.01451226
#> gs3 -0.19452894  0.21379668 -0.18038804 -0.13955592 -0.19972420
#> gs4  0.31206266 -0.02125529 -0.25247325  0.02560598 -0.33435424
#> gs5  0.18774756 -0.06777329 -0.02900350 -0.08483818  0.09954821
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
#>      gsva 2468.6860 2469.3363 2561.1002 2471.8374 2478.6492 3328.2424    10
#>  bixverse  480.4741  480.8475  481.8822  481.5707  483.1557  483.5452    10
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
#> gs1  0.13733024 -0.06231659  0.08237475  0.02371160  0.17354557
#> gs2  0.01141301 -0.11697659  0.15630981  0.03943545 -0.24538694
#> gs3 -0.01359398  0.17017295  0.19956091  0.06015560  0.17523035
#> gs4 -0.25416478  0.35893676  0.05329579 -0.06577850 -0.20355401
#> gs5  0.27387664  0.03733247 -0.15677380  0.22644909  0.07011793
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
#> gs1  0.13733024 -0.06231659  0.08237475  0.02371160  0.17354557
#> gs2  0.01141301 -0.11697659  0.15630981  0.03943545 -0.24538694
#> gs3 -0.01359398  0.17017295  0.19956091  0.06015560  0.17523035
#> gs4 -0.25416478  0.35893676  0.05329579 -0.06577850 -0.20355401
#> gs5  0.27387664  0.03733247 -0.15677380  0.22644909  0.07011793
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
#>      gsva 18.841476 18.846712 18.855018 18.851386 18.854761 18.887556    10
#>  bixverse  2.689627  2.692918  2.698892  2.693782  2.697896  2.732564    10
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
#>             s1         s2          s3          s4          s5
#> gs1 0.14739663 0.17154602  0.09733903  0.02210756  0.14363567
#> gs2 0.14121745 0.11569797  0.05234644 -0.02224332  0.15213105
#> gs3 0.02883741 0.27685852  0.02157257  0.02381283  0.08503531
#> gs4 0.27242318 0.13781409 -0.07864781  0.11931302 -0.06123245
#> gs5 0.24280877 0.07610833  0.05771482  0.08779286  0.16066583
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
#>             s1         s2          s3          s4          s5
#> gs1 0.14739663 0.17154602  0.09733903  0.02210756  0.14363567
#> gs2 0.14121745 0.11569797  0.05234644 -0.02224332  0.15213105
#> gs3 0.02883741 0.27685852  0.02157257  0.02381283  0.08503531
#> gs4 0.27242318 0.13781409 -0.07864781  0.11931302 -0.06123245
#> gs5 0.24280877 0.07610833  0.05771482  0.08779286  0.16066583
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
#>      gsva 620.6799 626.9014 720.8176 636.0647 654.0252 1484.4143    10
#>  bixverse 107.4565 107.7676 109.1566 108.3362 111.4066  111.9369    10
```
