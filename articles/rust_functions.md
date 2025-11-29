# Bixverse R \<\> Rust interface

## Using the Rust functions in more general R code

`bixverse` has a number of R \<\> Rust functions that might be useful
more generally. This vignette is there to explore how to use these
functions more general. This vignette will show you to use:

- Correlations (via [`rs_cor()`](../reference/rs_cor.md)).
- Covariance (via [`rs_covariance()`](../reference/rs_covariance.md)).
- Covariance to correlation (via
  [`rs_cov2cor()`](../reference/rs_cov2cor.md)).
- Distance calculations (vs [`rs_dist()`](../reference/rs_dist.md)).
- Correlations between two matrices (via
  [`rs_cor2()`](../reference/rs_cor2.md)).
- Calculate the mutual information between columns (via
  `rs_mutual_information()`)
- Calculations of set similarities between (via
  [`rs_set_similarity_list()`](../reference/rs_set_similarity_list.md)
  or [`rs_set_similarity()`](../reference/rs_set_similarity.md)).

``` r
library(bixverse)
```

This vignette will explore the interfaces in `bixverse` to underlying
Rust code you can leverage to make your code faster. But be careful… Any
function with `rs_` can cause panics when misused.

### Correlations, co-variance, and cosine distance

A typical thing one uses a lot are correlations, co-variance and cosine
distances. Let’s look at how to use `bixverse` to calculate these.

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

Let’s look at Pearson’s correlations.

``` r
r_pearsons_res <- cor(random_data)

rust_pearsons_res <- rs_cor(random_data, spearman = FALSE)

# tiny numerical differences do exist
all.equal(r_pearsons_res, rust_pearsons_res, tolerance = 1e-15)
#> [1] TRUE
```

Speed differences:

``` r
microbenchmark::microbenchmark(
  r = cor(random_data),
  rust = rs_cor(random_data, spearman = FALSE),
  times = 50L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median        uq      max neval
#>     r 134.124470 134.267096 134.636403 134.339521 134.50529 140.8720    50
#>  rust   4.942054   5.609561   7.318338   7.160285   8.22372  15.1611    50
```

#### Spearman’s correlations

Let’s look at Spearman’s correlations…

``` r
r_spearman_res <- cor(random_data, method = "spearman")

rust_spearman_res <- rs_cor(random_data, spearman = TRUE)

# tiny numerical differences do exist
all.equal(r_spearman_res, rust_spearman_res, tolerance = 1e-15)
#> [1] TRUE
```

Speed differences:

``` r
microbenchmark::microbenchmark(
  r = cor(random_data, method = "spearman"),
  rust = rs_cor(random_data, spearman = TRUE),
  times = 50L
)
#> Unit: milliseconds
#>  expr        min        lq       mean     median        uq       max neval
#>     r 164.163660 167.14326 170.979253 167.657636 168.36063 250.75675    50
#>  rust   6.776269   7.16816   7.891583   7.445678   8.01459  13.47219    50
```

#### Covariance

And covariance calculations

``` r
r_covar_res <- cov(random_data)

rust_covar_res <- rs_covariance(random_data)

# tiny numerical differences do exist
all.equal(r_covar_res, rust_covar_res, tolerance = 1e-15)
#> [1] TRUE
```

Speed differences:

``` r
microbenchmark::microbenchmark(
  r = cov(random_data),
  rust = rs_covariance(random_data),
  times = 50L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median        uq       max neval
#>     r 133.857461 133.950575 134.384759 134.015792 134.13541 138.35896    50
#>  rust   4.179711   4.424809   5.615987   5.173211   5.92514  11.53455    50
```

Let’s transform covariance to correlations

``` r
r_cor_from_covar <- cov2cor(r_covar_res)

rust_cor_from_covar <- rs_cov2cor(rust_covar_res)

# tiny numerical differences do exist
all.equal(r_cor_from_covar, rust_cor_from_covar, tolerance = 1e-15)
#> [1] TRUE
```

How much faster is Rust here? (This is where we observe the least
improvement.)

``` r
microbenchmark::microbenchmark(
  r = cov2cor(r_covar_res),
  rust = rs_cov2cor(rust_covar_res),
  times = 50L
)
#> Unit: milliseconds
#>  expr      min       lq     mean   median       uq      max neval
#>     r 2.090963 2.144983 2.784570 2.167902 2.249339 6.970802    50
#>  rust 1.377772 1.410663 1.637771 1.470384 1.514106 5.424866    50
```

#### Correlations between two matrices

We also provide interface to do correlations between two sets of data

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

Let’s test the correlations between the two matrices

``` r
r_cor_2_matrices <- cor(random_data, random_data_2)

rust_cor_2_matrices <- rs_cor2(random_data, random_data_2, spearman = FALSE)

# The small precision differences are further propagated in here
all.equal(r_cor_2_matrices, rust_cor_2_matrices, tolerance = 1e-14)
#> [1] TRUE
```

And speed between Rust and R

``` r
microbenchmark::microbenchmark(
  r = cor(random_data, random_data_2),
  rust = rs_cor2(random_data, random_data_2, spearman = FALSE),
  times = 50L
)
#> Unit: milliseconds
#>  expr        min         lq     mean     median         uq       max neval
#>     r 190.622887 190.770412 190.9804 190.881664 190.984061 195.54400    50
#>  rust   4.743404   5.351489   6.2463   5.550019   6.958268  13.60956    50
```

### Distance metrics

`bixverse` also includes various distance metrics. Let’s have a look at
very common ones.

#### Euclidean distance

Let’s check out the very classical Euclidean distance.

``` r
r_euclidean_distance <- as.matrix(
  dist(
    random_data, 
    method = "euclidean"
  )
)

rs_euclidean_distance <- rs_dist(
  t(random_data), # dist() calculates row-wise distances
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

Again, Rust is much faster…

``` r
# to remove any overhead, we do not do the matrix transformation
# or the transpositation. due to the equal dimensions, it will
# not matter
microbenchmark::microbenchmark(
  r = dist(random_data, method = "euclidean"),
  rust = rs_dist(random_data, distance_type = "euclidean"),
  times = 50L
)
#> Unit: milliseconds
#>  expr       min       lq      mean    median        uq      max neval
#>     r 96.404274 96.64864 96.792907 96.766674 96.925937 97.30067    50
#>  rust  4.590238  5.34705  6.108914  5.610067  6.874722 10.97611    50
```

#### Other distance metrics

Let’s test other distance metrics…

**Manhattan distance**

We observe again equivalence of Rust and R results.

``` r
r_manhattan_distance <- as.matrix(
  dist(
    random_data, 
    method = "manhattan"
  )
)

rs_manhattan_distance <- rs_dist(
  t(random_data), # dist() calculates row-wise distances
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

And speed ups in Rust.

``` r
# to remove any overhead, we do not do the matrix transformation
# or the transpositation. due to the equal dimensions, it will
# not matter
microbenchmark::microbenchmark(
  r = dist(random_data, method = "manhattan"),
  rust = rs_dist(random_data, distance_type = "manhattan"),
  times = 50L
)
#> Unit: milliseconds
#>  expr      min       lq     mean   median       uq      max neval
#>     r 90.91773 91.10968 91.44102 91.28737 91.42767 96.64115    50
#>  rust 17.05988 17.30768 17.72751 17.57983 17.92434 21.21891    50
```

**Canberra distance**

Equivalence observed…

``` r
r_canberra_distance <- as.matrix(
  dist(
    random_data, 
    method = "canberra"
  )
)

rs_canberra_distance <- rs_dist(
  t(random_data), # dist() calculates row-wise distances
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

… and speed-ups.

``` r
# to remove any overhead, we do not do the matrix transformation
# or the transpositation. due to the equal dimensions, it will
# not matter
microbenchmark::microbenchmark(
  r = dist(random_data, method = "canberra"),
  rust = rs_dist(random_data, distance_type = "canberra"),
  times = 50L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median        uq       max neval
#>     r 136.69158 137.13483 137.41841 137.29600 137.47533 139.98934    50
#>  rust  45.81207  46.01443  46.19196  46.19184  46.35743  46.62335    50
```

### Mutual information

Another quantity you might be interested in is the mutual information
between between your variables of interest. `bixverse` also provides a
Rust-powered interface into very rapid calculations in this regard.

``` r
# we are going to create a reduced set for the benchmark to finish 
# in a reasonable time
set.seed(246L)

nrows <- 100
ncols <- 100

mat <- matrix(data = rnorm(nrows * ncols), nrow = nrows, ncol = ncols)
rownames(mat) <- sprintf("sample_%i", 1:nrows)
colnames(mat) <- sprintf("feature_%i", 1:ncols)
```

To run the function you can just use the following code. The Rust
version will default to sqrt(nrow())

``` r
rust_res_mi <- rs_mutual_info(
  mat, 
  n_bins = NULL, 
  normalise = FALSE,
  strategy = "equal_width"
)
rownames(rust_res_mi) <- colnames(rust_res_mi) <- colnames(mat)

# ensure that the same discretisation is used
infotheo_res_mi <- infotheo::mutinformation(infotheo::discretize(
  mat,
  disc = "equalwidth",
  nbins = sqrt(nrow(mat))
))

all.equal(rust_res_mi, infotheo_res_mi)
#> [1] TRUE
```

There is also an equal frequency strategy implemented:

``` r
rust_res_mi <- rs_mutual_info(
  mat, 
  n_bins = NULL, 
  normalise = FALSE,
  strategy = "equal_freq"
)
rownames(rust_res_mi) <- colnames(rust_res_mi) <- colnames(mat)

# ensure that the same discretisation is used
infotheo_res_mi <- infotheo::mutinformation(infotheo::discretize(
  mat,
  disc = "equalfreq",
  nbins = sqrt(nrow(mat))
))

all.equal(rust_res_mi, infotheo_res_mi)
#> [1] TRUE
```

Again, the Rust implementation is way faster.

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
  times = 50L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq       max neval
#>  infotheo 82.751570 83.530403 84.201699 84.084402 84.533626 89.944796    50
#>      rust  3.354581  3.405336  3.581643  3.466249  3.767241  3.927892    50
```

And with equal frequency per bin:

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
  times = 50L
)
#> Unit: milliseconds
#>      expr        min         lq       mean     median         uq        max
#>  infotheo 100.368814 100.798355 101.550881 101.178675 101.715817 114.053637
#>      rust   4.383572   4.429909   4.600883   4.577943   4.722164   5.108525
#>  neval
#>     50
#>     50
```

### Set similarities

`bixverse` also provides various functions that leverage Rust to do set
similarities. Let’s compare this against some R implementation. Below is
a naive R implementation to calculate the Jaccard similarities between
two lists containing strings.

``` r
jaccard_sim <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

set.seed(42L)

random_sets_1 <- purrr::map(1:1000L, ~{
  paste(sample(LETTERS, sample(1:20, 1)))
})

random_sets_2 <- purrr::map(1:500L, ~{
  paste(sample(LETTERS, sample(1:20, 1)))
})
```

Let’s start with a naive [purrr](https://purrr.tidyverse.org) approach.
The problem obviously will benefit from some form of parallel
processing, but for completeness, let’s run this version, too.

``` r
tictoc::tic()

r_results <- purrr::map(random_sets_1, \(x) {
  purrr::map_dbl(random_sets_2, \(y) {jaccard_sim(x, y)})
}, .progress = TRUE)
#>  ■■■■                              11% |  ETA: 23s
#>  ■■■■■■■■                          23% |  ETA: 20s
#>  ■■■■■■■■■■■                       35% |  ETA: 16s
#>  ■■■■■■■■■■■■■■■                   46% |  ETA: 14s
#>  ■■■■■■■■■■■■■■■■■■■               58% |  ETA: 10s
#>  ■■■■■■■■■■■■■■■■■■■■■■            70% |  ETA:  7s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■        82% |  ETA:  5s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% |  ETA:  2s

similarity_matrix <- matrix(
  data = unlist(r_results), 
  nrow = length(random_sets_1), 
  ncol = length(random_sets_2),
  byrow = TRUE
)
tictoc::toc()
#> 25.244 sec elapsed
```

Let’s try to parallelise this vi [furrr](https://furrr.futureverse.org)
and evaluate how much we can improve this.

``` r
future::plan(strategy = future::multisession(workers = parallel::detectCores()))

tictoc::tic()

r_results <- furrr::future_map(random_sets_1, \(x) {
  purrr::map_dbl(random_sets_2, \(y) {jaccard_sim(x, y)})
}, .progress = TRUE)

similarity_matrix <- matrix(
  data = unlist(r_results), 
  nrow = length(random_sets_1), 
  ncol = length(random_sets_2),
  byrow = TRUE
)
tictoc::toc()
#> 13.99 sec elapsed

future::plan(strategy = future::sequential())
```

Parallelisation as expected accelerates the calculations substantially.
What about [mirai](https://mirai.r-lib.org/articles/mirai.html)? It
should produce less overhead compared to future.

``` r
mirai::daemons(parallel::detectCores())

tictoc::tic()

r_results_mirai <- mirai::mirai_map(random_sets_1, \(x, sets_2) {
  purrr::map_dbl(sets_2, \(y) {
    length(intersect(x, y)) / length(union(x, y))
  })
}, .args = list(sets_2 = random_sets_2))[]

similarity_matrix <- matrix(
  data = unlist(r_results_mirai), 
  nrow = length(random_sets_1), 
  ncol = length(random_sets_2),
  byrow = TRUE
)

tictoc::toc()
#> 13.315 sec elapsed

mirai::daemons(0)
```

That’s already a bit better. Maybe we can optimise the Jaccard
similarity function further and accelerate this even more.

``` r
tictoc::tic()
mirai::daemons(4)

r_results_mirai <- mirai::mirai_map(random_sets_1, \(x, sets_2) {
  purrr::map_dbl(sets_2, \(y) {
    # Optimised function
    intersection <- length(intersect(x, y))
    union_size <- length(x) + length(y) - intersection
    intersection / union_size
  })
}, .args = list(sets_2 = random_sets_2))[]

similarity_matrix <- matrix(
  data = unlist(r_results), 
  nrow = length(random_sets_1), 
  ncol = length(random_sets_2),
  byrow = TRUE
)

mirai::daemons(0)
tictoc::toc()
#> 8.081 sec elapsed
```

This is already much better… Pending on your system you might be seeing
a 4x increase in speed compared to the naive purrr version. However,
let’s compare against the Rust version implemented in `bixverse`:

``` r
tictoc::tic()
rust_res <- rs_set_similarity_list2(
    s_1_list = random_sets_1,
    s_2_list = random_sets_2,
  overlap_coefficient = FALSE
)
tictoc::toc()
#> 0.049 sec elapsed

all.equal(similarity_matrix, rust_res, tolerance = 1e-15)
#> [1] TRUE
```

Yep, that is MUCH faster… These are some of the functions that are
exposed in `bixverse` and can be integrated into other workflows. There
are many more highly specialised functions that can be used in other
workflows. It’s best to explore the package and the underlying code base
to identify all of these functions.
