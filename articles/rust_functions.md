# Bixverse R \<\> Rust interface

## Why Rust?

Honest answer? The creation of this package was based on the desire of
the author to play with Rust. To make it useful, the initial playground
was stuff relevant to day-to-day work in bioinformatics and
computational biology: correlations, PCA, other matrix factorisation
methods, and so on. It quickly became apparent that Rust delivered
incredible speed gains over R — sometimes 100x — with millions of
hypergeometric tests completing in mere seconds. The results were
striking. Memory safety and the fantastic work by the [rextendr
team](https://github.com/extendr/rextendr) made interfacing R and Rust
straightforward and, frankly, quite fun.

### Correlations and covariance

One of the most typical operations in computational biology is computing
correlations and covariances; and these were the first functions
implemented, largely to understand
[faer](https://github.com/sarah-ek/faer-rs) and the Rust/R interface.
What becomes apparent very quickly is the following. For 1000 x 1000
matrices on an M1 Max MacBook Pro (using Apple’s Accelerate framework
for R BLAS):

**Covariance**

| Implementation | Median (ms, 20 runs) |
|---------------:|---------------------:|
|         base R |               515.73 |
|           Rust |                18.66 |

*A 27x increase in speed.*

**Pearson correlation**

| Implementation | Median (ms, 20 runs) |
|---------------:|---------------------:|
|         base R |               512.83 |
|           Rust |                20.90 |

*A 24x increase in speed.*

**Spearman correlation**

| Implementation | Median (ms, 20 runs) |
|---------------:|---------------------:|
|         base R |               603.00 |
|           Rust |                29.47 |

*A 20x increase in speed.*

The gains are even more pronounced on larger matrices. For a 1000 x 5000
matrix:

**Pearson correlation (1000 x 5000)**

| Implementation | Median (ms, 20 runs) |
|---------------:|---------------------:|
|         base R |             12994.37 |
|           Rust |               419.88 |

*A 30x increase in speed.*

Given that matrix algebra underlies a large proportion of computational
bioinformatics workflows, this makes a compelling case for Rust under
the hood. The rest of this vignette shows how to use the exposed Rust
functions directly, with runnable benchmarks you can verify on your own
machine.

## Using the Rust functions in more general R code

`bixverse` exposes a number of R/Rust functions that are useful beyond
the higher-level package functionality. This vignette covers:

- Correlations (via
  [`rs_cor()`](https://gregorlueg.github.io/bixverse/reference/rs_cor.md)).
- Covariance (via
  [`rs_covariance()`](https://gregorlueg.github.io/bixverse/reference/rs_covariance.md)).
- Covariance to correlation (via
  [`rs_cov2cor()`](https://gregorlueg.github.io/bixverse/reference/rs_cov2cor.md)).
- Distance calculations (via
  [`rs_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_dist.md)).
- Correlations between two matrices (via
  [`rs_cor2()`](https://gregorlueg.github.io/bixverse/reference/rs_cor2.md)).
- Mutual information between columns (via `rs_mutual_information()`).
- Set similarities (via
  [`rs_set_similarity_list()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity_list.md)
  or
  [`rs_set_similarity()`](https://gregorlueg.github.io/bixverse/reference/rs_set_similarity.md)).

These may be useful for your own package development.

> **Note**
>
> The vignette was built on a GitHub runner which is a 2-core machine
> with ~8 GB of memory. The performance differences tend to be (even)
> more pronounced on machines with more cores and computational power.

## Set-up

``` r

library(bixverse)
```

Any function prefixed with `rs_` calls directly into Rust and can panic
if misused — check inputs before passing them through. There is a LOT of
these functions in the package, and they are all exported. Why? I want
other people to be able to integrate the code into their own R packages,
or if they feel brave, they can just use the Rust crate directly, see
[here](https://crates.io/crates/bixverse-rs).

### Correlations, covariance, and cosine distance

``` r

no_rows <- 500L
no_cols <- 500L

set.seed(10101L)

random_data <- matrix(
  data = rnorm(no_rows * no_cols),
  ncol = no_cols,
  nrow = no_rows
)
```

#### Pearson’s correlation

``` r

r_pearsons_res <- cor(random_data)

rust_pearsons_res <- rs_cor(random_data, spearman = FALSE)

# tiny numerical differences do exist
all.equal(r_pearsons_res, rust_pearsons_res, tolerance = 1e-15)
#> [1] TRUE
```

Let’s check the speed.

``` r

microbenchmark::microbenchmark(
  r = cor(random_data),
  rust = rs_cor(random_data, spearman = FALSE),
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean    median         uq       max neval
#>     r 110.031310 110.129136 110.354116 110.29647 110.452717 111.16280    10
#>  rust   2.743614   2.771826   3.600472   3.67507   4.191308   4.41399    10
```

#### Spearman’s correlation

What about ranked correlations?

``` r

r_spearman_res <- cor(random_data, method = "spearman")

rust_spearman_res <- rs_cor(random_data, spearman = TRUE)

all.equal(r_spearman_res, rust_spearman_res, tolerance = 1e-15)
#> [1] TRUE
```

And speed

``` r

microbenchmark::microbenchmark(
  r = cor(random_data, method = "spearman"),
  rust = rs_cor(random_data, spearman = TRUE),
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq        max neval
#>     r 135.096778 136.514829 137.905242 138.546487 139.173305 139.736934    10
#>  rust   4.199831   4.881545   5.023703   5.029326   5.115835   6.499484    10
```

#### Covariance

Covariance matrices are also very common…

``` r

r_covar_res <- cov(random_data)

rust_covar_res <- rs_covariance(random_data)

all.equal(r_covar_res, rust_covar_res, tolerance = 1e-15)
#> [1] TRUE
```

And speed comparisons

``` r

microbenchmark::microbenchmark(
  r = cov(random_data),
  rust = rs_covariance(random_data),
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean    median         uq        max neval
#>     r 109.784222 109.832164 109.865821 109.84001 109.867337 110.103347    10
#>  rust   1.729474   1.822843   1.864386   1.86126   1.937945   1.987308    10
```

#### Covariance to correlation

``` r

r_cor_from_covar <- cov2cor(r_covar_res)

rust_cor_from_covar <- rs_cov2cor(rust_covar_res)

all.equal(r_cor_from_covar, rust_cor_from_covar, tolerance = 1e-15)
#> [1] TRUE
```

This is where the Rust implementation shows the least improvement — the
operation is already cheap in base R — but the interface remains
consistent with the rest of the `rs_` family.

``` r

microbenchmark::microbenchmark(
  r = cov2cor(r_covar_res),
  rust = rs_cov2cor(rust_covar_res),
  times = 10L
)
#> Unit: microseconds
#>  expr      min       lq      mean   median       uq      max neval
#>     r 1778.197 1787.531 2083.1887 1808.812 1822.903 4527.398    10
#>  rust  207.909  220.789  605.1238  236.597  268.409 3776.441    10
```

#### Correlations between two matrices

So far, all of the functions did pairwise column-based correlations.
What about two correlations against each other?

``` r

no_rows_2 <- 500L
no_cols_2 <- 400L

set.seed(23L)

random_data_2 <- matrix(
  data = rnorm(no_rows_2 * no_cols_2, mean = 2, sd = 2),
  ncol = no_cols_2,
  nrow = no_rows_2
)
```

Let’s compare against R:

``` r

r_cor_2_matrices <- cor(random_data, random_data_2)

rust_cor_2_matrices <- rs_cor2(random_data, random_data_2, spearman = FALSE)

# small precision differences are further propagated here
all.equal(r_cor_2_matrices, rust_cor_2_matrices, tolerance = 1e-14)
#> [1] TRUE
```

And speed:

``` r

microbenchmark::microbenchmark(
  r = cor(random_data, random_data_2),
  rust = rs_cor2(random_data, random_data_2, spearman = FALSE),
  times = 10L
)
#> Unit: milliseconds
#>  expr       min         lq       mean     median         uq        max neval
#>     r 180.94765 181.077626 181.137108 181.082588 181.279746 181.372625    10
#>  rust   2.97001   3.056008   3.227348   3.181741   3.273272   3.935368    10
```

### Distance metrics

#### Euclidean distance

Note that [`dist()`](https://rdrr.io/r/stats/dist.html) computes
row-wise distances, while
[`rs_dist()`](https://gregorlueg.github.io/bixverse/reference/rs_dist.md)
expects columns to represent observations — hence the transpose.

``` r

r_euclidean_distance <- as.matrix(
  dist(random_data, method = "euclidean")
)

rs_euclidean_distance <- rs_dist(
  t(random_data),
  distance_type = "euclidean"
)

all.equal(
  r_euclidean_distance,
  rs_euclidean_distance,
  tolerance = 1e-14,
  check.attributes = FALSE
)
#> [1] TRUE
```

And speed…

``` r

microbenchmark::microbenchmark(
  r = dist(random_data, method = "euclidean"),
  rust = rs_dist(random_data, distance_type = "euclidean"),
  times = 10L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median       uq       max neval
#>     r 82.446143 82.552101 82.736633 82.571685 82.86888 83.578669    10
#>  rust  1.799238  1.815572  1.890168  1.890529  1.95503  2.009611    10
```

#### Manhattan distance

``` r

r_manhattan_distance <- as.matrix(
  dist(random_data, method = "manhattan")
)

rs_manhattan_distance <- rs_dist(
  t(random_data),
  distance_type = "manhattan"
)

all.equal(
  r_manhattan_distance,
  rs_manhattan_distance,
  tolerance = 1e-14,
  check.attributes = FALSE
)
#> [1] TRUE
```

And speed… ?

``` r

microbenchmark::microbenchmark(
  r = dist(random_data, method = "manhattan"),
  rust = rs_dist(random_data, distance_type = "manhattan"),
  times = 10L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median        uq       max neval
#>     r 73.740455 73.955095 74.034302 74.011344 74.071848 74.421249    10
#>  rust  4.879412  4.924549  5.118144  5.037683  5.290765  5.603841    10
```

#### Canberra distance

``` r

r_canberra_distance <- as.matrix(
  dist(random_data, method = "canberra")
)

rs_canberra_distance <- rs_dist(
  t(random_data),
  distance_type = "canberra"
)

all.equal(
  r_canberra_distance,
  rs_canberra_distance,
  tolerance = 1e-14,
  check.attributes = FALSE
)
#> [1] TRUE
```

Same benchmarks as before…

``` r

microbenchmark::microbenchmark(
  r = dist(random_data, method = "canberra"),
  rust = rs_dist(random_data, distance_type = "canberra"),
  times = 10L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median        uq       max neval
#>     r 110.86916 111.17974 111.64728 111.27110 111.57319 113.99811    10
#>  rust  11.02146  11.06175  11.21033  11.19129  11.27857  11.62513    10
```

### Mutual information

``` r

set.seed(246L)

nrows <- 100
ncols <- 100

mat <- matrix(data = rnorm(nrows * ncols), nrow = nrows, ncol = ncols)
rownames(mat) <- sprintf("sample_%i", 1:nrows)
colnames(mat) <- sprintf("feature_%i", 1:ncols)
```

The Rust version defaults to `sqrt(nrow())` bins when `n_bins = NULL`.
Two discretisation strategies are available.

**Equal width:**

``` r

rust_res_mi <- rs_mutual_info(
  mat,
  n_bins = NULL,
  normalise = FALSE,
  strategy = "equal_width"
)
rownames(rust_res_mi) <- colnames(rust_res_mi) <- colnames(mat)

infotheo_res_mi <- infotheo::mutinformation(infotheo::discretize(
  mat,
  disc = "equalwidth",
  nbins = sqrt(nrow(mat))
))

all.equal(rust_res_mi, infotheo_res_mi)
#> [1] TRUE
```

``` r

microbenchmark::microbenchmark(
  infotheo = infotheo::mutinformation(infotheo::discretize(
    mat,
    disc = "equalwidth",
    nbins = sqrt(nrow(mat))
  )),
  rust = rs_mutual_info(
    mat,
    n_bins = NULL,
    normalise = FALSE,
    strategy = "equal_width"
  ),
  times = 10L
)
#> Unit: milliseconds
#>      expr      min       lq      mean    median        uq       max neval
#>  infotheo 74.19109 74.80354 75.250287 75.199142 75.720292 76.734782    10
#>      rust  1.92805  2.04925  2.156806  2.101943  2.170391  2.666078    10
```

**Equal frequency:**

``` r

rust_res_mi <- rs_mutual_info(
  mat,
  n_bins = NULL,
  normalise = FALSE,
  strategy = "equal_freq"
)
rownames(rust_res_mi) <- colnames(rust_res_mi) <- colnames(mat)

infotheo_res_mi <- infotheo::mutinformation(infotheo::discretize(
  mat,
  disc = "equalfreq",
  nbins = sqrt(nrow(mat))
))

all.equal(rust_res_mi, infotheo_res_mi)
#> [1] TRUE
```

``` r

microbenchmark::microbenchmark(
  infotheo = infotheo::mutinformation(infotheo::discretize(
    mat,
    disc = "equalfreq",
    nbins = sqrt(nrow(mat))
  )),
  rust = rs_mutual_info(
    mat,
    n_bins = NULL,
    normalise = FALSE,
    strategy = "equal_freq"
  ),
  times = 10L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq       max neval
#>  infotheo 91.409173 92.383904 93.446272 92.930773 93.967220 96.790308    10
#>      rust  2.683565  2.698597  2.860279  2.811445  2.848209  3.659527    10
```

### Set similarities

`bixverse` also exposes Rust-accelerated set similarity calculations. To
illustrate the performance, below is a comparison against a naive R
implementation of Jaccard similarity across two lists.

``` r

jaccard_sim <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

set.seed(42L)

random_sets_1 <- purrr::map(
  1:1000L,
  ~ paste(sample(LETTERS, sample(1:20, 1)))
)

random_sets_2 <- purrr::map(
  1:500L,
  ~ paste(sample(LETTERS, sample(1:20, 1)))
)
```

Starting with a sequential [purrr](https://purrr.tidyverse.org) approach
as the baseline:

``` r

tictoc::tic()

r_results <- purrr::map(
  random_sets_1,
  \(x) {
    purrr::map_dbl(random_sets_2, \(y) jaccard_sim(x, y))
  },
  .progress = TRUE
)
#>  ■■■■■                             13% |  ETA: 18s
#>  ■■■■■■■■■■                        29% |  ETA: 14s
#>  ■■■■■■■■■■■■■■                    44% |  ETA: 11s
#>  ■■■■■■■■■■■■■■■■■■■               60% |  ETA:  8s
#>  ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  5s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA:  2s

similarity_matrix <- matrix(
  data = unlist(r_results),
  nrow = length(random_sets_1),
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 19.718 sec elapsed
```

Parallelising via [furrr](https://furrr.futureverse.org) gives a
meaningful speedup:

``` r

future::plan(strategy = future::multisession(workers = parallel::detectCores()))

tictoc::tic()

r_results <- furrr::future_map(
  random_sets_1,
  \(x) {
    purrr::map_dbl(random_sets_2, \(y) jaccard_sim(x, y))
  },
  .progress = TRUE
)

similarity_matrix <- matrix(
  data = unlist(r_results),
  nrow = length(random_sets_1),
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 9.377 sec elapsed

future::plan(strategy = future::sequential())
```

[mirai](https://mirai.r-lib.org/articles/mirai.html) reduces process
scheduling overhead (somtimes) compared to `future` and does a bit
better:

``` r

mirai::daemons(parallel::detectCores())

tictoc::tic()

r_results_mirai <- mirai::mirai_map(
  random_sets_1,
  \(x, sets_2) {
    purrr::map_dbl(sets_2, \(y) {
      length(intersect(x, y)) / length(union(x, y))
    })
  },
  .args = list(sets_2 = random_sets_2)
)[]

similarity_matrix <- matrix(
  data = unlist(r_results_mirai),
  nrow = length(random_sets_1),
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 9.489 sec elapsed

mirai::daemons(0)
```

Depending on your system this optimised parallel R version may be around
4x faster than the naive sequential baseline. Now compare against the
Rust implementation:

``` r

tictoc::tic()

rust_res <- rs_set_similarity_list2(
  s_1_list = random_sets_1,
  s_2_list = random_sets_2,
  overlap_coefficient = FALSE
)

tictoc::toc()
#> 0.041 sec elapsed

all.equal(similarity_matrix, rust_res, tolerance = 1e-15)
#> [1] TRUE
```

The margin here is not subtle. The functions shown in this vignette are
a representative sample of what is exposed in `bixverse`: there are many
more specialised functions throughout the package worth exploring (if
you are brave and can deal with a panic here and there).
