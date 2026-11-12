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
#>  expr        min         lq      mean     median         uq        max neval
#>     r 110.062367 110.136107 110.32863 110.162760 110.370584 111.349815    10
#>  rust   2.783052   2.892474   3.64625   3.741562   3.920978   4.480046    10
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
#>     r 135.022567 136.991453 138.501960 138.846356 139.241746 143.263843    10
#>  rust   4.302443   4.370653   5.140332   5.090354   5.370064   6.580978    10
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
#>  expr        min         lq       mean     median         uq        max neval
#>     r 109.780240 109.795983 109.910812 109.878505 110.035738 110.084120    10
#>  rust   1.784313   1.840496   1.908916   1.910119   1.971831   2.028194    10
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
#>  expr      min       lq      mean    median       uq      max neval
#>     r 1776.552 1791.834 2083.0815 1817.0915 1828.819 4492.274    10
#>  rust  193.637  224.032  511.6018  238.0225  251.883 2981.086    10
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
#>  expr        min         lq       mean     median         uq        max neval
#>     r 181.059704 181.102948 181.646932 181.342813 181.891326 183.704952    10
#>  rust   2.937762   3.030349   3.227149   3.203761   3.307038   3.811246    10
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
#>  expr       min        lq      mean    median        uq       max neval
#>     r 82.815040 82.993926 83.245662 83.115375 83.378485 84.142376    10
#>  rust  1.825073  1.860386  1.956383  1.957289  2.051809  2.074052    10
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
#>     r 73.949219 74.020564 74.570372 74.483084 74.710496 76.205711    10
#>  rust  4.883553  5.032364  5.211405  5.096544  5.180093  6.071783    10
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
#>     r 108.91149 110.78884 110.98411 111.04993 111.32994 112.25427    10
#>  rust  11.05797  11.10622  11.19424  11.17402  11.28564  11.36503    10
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
#>  infotheo 76.217699 76.421491 76.860125 76.730218 77.002281 77.982924    10
#>      rust  1.911021  1.935607  2.086539  2.041543  2.160039  2.630146    10
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
#>      expr       min        lq      mean   median        uq       max neval
#>  infotheo 94.702407 94.896214 95.331694 95.47893 95.731772 95.943115    10
#>      rust  2.663035  2.686188  2.766655  2.80164  2.829491  2.844804    10
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
#>  ■■■■■                             14% |  ETA: 18s
#>  ■■■■■■■■■                         28% |  ETA: 15s
#>  ■■■■■■■■■■■■■■                    44% |  ETA: 12s
#>  ■■■■■■■■■■■■■■■■■■■               59% |  ETA:  8s
#>  ■■■■■■■■■■■■■■■■■■■■■■■           74% |  ETA:  5s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% |  ETA:  2s

similarity_matrix <- matrix(
  data = unlist(r_results),
  nrow = length(random_sets_1),
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 20.164 sec elapsed
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
#> 9.352 sec elapsed

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
#> 9.52 sec elapsed

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
#> 0.044 sec elapsed

all.equal(similarity_matrix, rust_res, tolerance = 1e-15)
#> [1] TRUE
```

The margin here is not subtle. The functions shown in this vignette are
a representative sample of what is exposed in `bixverse`: there are many
more specialised functions throughout the package worth exploring (if
you are brave and can deal with a panic here and there).
