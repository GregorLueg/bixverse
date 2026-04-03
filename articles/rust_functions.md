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
#>  expr        min         lq       mean     median         uq        max neval
#>     r 135.312394 135.373840 135.431313 135.439241 135.495697 135.531073    10
#>  rust   5.444358   5.592837   6.286816   5.783357   6.710764   8.378047    10
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
#>  expr        min         lq       mean     median         uq       max neval
#>     r 165.588371 166.797249 168.556694 168.947550 170.386463 170.90678    10
#>  rust   6.726233   7.417163   8.263062   7.587611   9.250427  10.84584    10
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
#>  expr       min         lq     mean   median         uq        max neval
#>     r 134.96988 135.093967 135.1526 135.1751 135.198401 135.252323    10
#>  rust   4.79169   4.952139   5.6938   5.0513   5.679457   8.769189    10
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
#> Unit: milliseconds
#>  expr      min       lq     mean   median       uq      max neval
#>     r 2.207443 2.357042 3.428634 3.269662 3.914912 6.063334    10
#>  rust 1.425653 1.595239 2.871570 2.070142 3.600184 7.167337    10
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
#>  expr        min         lq       mean     median         uq      max neval
#>     r 193.945214 193.992753 194.097937 194.052764 194.209998 194.3461    10
#>  rust   4.752236   5.378115   6.090929   5.801671   6.629693   7.8007    10
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
#>  expr      min       lq      mean    median        uq       max neval
#>     r 97.05414 97.12536 97.454028 97.193399 97.286984 99.144955    10
#>  rust  4.55638  4.65266  4.858743  4.849523  5.068346  5.141383    10
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
#>  expr      min       lq     mean   median       uq      max neval
#>     r 91.67173 91.74867 91.93962 91.93518 92.11489 92.33481    10
#>  rust 17.18661 17.34389 17.51150 17.45779 17.63719 18.01073    10
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
#>  expr       min       lq     mean    median        uq       max neval
#>     r 137.55280 137.5934 137.9984 137.72204 138.01359 140.17053    10
#>  rust  46.05706  46.2259  47.1314  46.58501  46.69156  51.13179    10
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
#>      expr       min        lq      mean    median        uq       max neval
#>  infotheo 84.479393 85.042795 85.288859 85.201020 85.589156 86.101162    10
#>      rust  3.324189  3.387938  3.519017  3.550006  3.654415  3.674012    10
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
#>      expr        min         lq       mean     median         uq        max
#>  infotheo 102.158603 102.618111 102.979135 103.064695 103.303240 103.561143
#>      rust   4.347139   4.430725   4.585577   4.607079   4.735558   4.777984
#>  neval
#>     10
#>     10
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
#>  ■■■                                7% |  ETA: 23s
#>  ■■■■■■■                           18% |  ETA: 20s
#>  ■■■■■■■■■■                        31% |  ETA: 17s
#>  ■■■■■■■■■■■■■■                    43% |  ETA: 14s
#>  ■■■■■■■■■■■■■■■■■■                55% |  ETA: 11s
#>  ■■■■■■■■■■■■■■■■■■■■■             67% |  ETA:  8s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■         80% |  ETA:  5s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      92% |  ETA:  2s

similarity_matrix <- matrix(
  data = unlist(r_results),
  nrow = length(random_sets_1),
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 24.684 sec elapsed
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
#> 14.028 sec elapsed

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
#> 13.252 sec elapsed

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
#> 0.047 sec elapsed

all.equal(similarity_matrix, rust_res, tolerance = 1e-15)
#> [1] TRUE
```

The margin here is not subtle. The functions shown in this vignette are
a representative sample of what is exposed in `bixverse`: there are many
more specialised functions throughout the package worth exploring (if
you are brave and can deal with a panic here and there).
