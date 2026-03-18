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

They might be useful for your down package development.

## Set-up

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
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq        max neval
#>     r 135.116561 135.374703 135.694504 135.451325 135.615893 137.517553    10
#>  rust   4.833584   5.687348   6.611594   6.064698   7.519058   9.437831    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq       max neval
#>     r 165.934707 166.705857 168.693532 169.082249 170.122426 170.44173    10
#>  rust   6.830593   7.456691   8.556271   7.787109   9.668221  11.24746    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq        max neval
#>     r 134.807513 134.909675 135.723772 135.164806 137.266274 137.543692    10
#>  rust   4.410174   5.248078   5.956568   5.339143   6.164369   9.031191    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr      min       lq     mean   median       uq      max neval
#>     r 2.202402 2.600345 3.587493 3.106631 4.725663 6.442809    10
#>  rust 1.416324 1.621778 3.063695 2.211825 4.044251 7.315849    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq        max neval
#>     r 193.675318 193.788000 193.943031 193.861592 194.090434 194.336814    10
#>  rust   4.885381   5.664866   6.639908   6.452732   7.354741   8.657143    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median       uq       max neval
#>     r 96.781609 97.136392 97.604174 97.455328 97.81650 99.277270    10
#>  rust  4.661733  4.693232  5.096851  5.049387  5.52716  5.614954    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr      min       lq     mean   median       uq      max neval
#>     r 91.70112 91.82370 92.22460 92.10353 92.67924 92.98861    10
#>  rust 17.37504 17.71626 17.88436 17.89619 18.14033 18.27615    10
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
  times = 10L
)
#> Unit: milliseconds
#>  expr       min        lq      mean    median        uq       max neval
#>     r 137.31488 137.44399 137.77648 137.60419 137.89925 139.28767    10
#>  rust  46.04405  46.33603  46.46993  46.46933  46.66423  46.80639    10
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
  times = 10L
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq     max neval
#>  infotheo 83.936368 84.467209 84.731817 84.861846 85.008380 85.1825    10
#>      rust  3.355194  3.383918  3.623576  3.669992  3.755551  4.0067    10
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
  times = 10L
)
#> Unit: milliseconds
#>      expr        min        lq       mean     median         uq        max
#>  infotheo 101.078833 101.28361 101.527034 101.530731 101.755696 102.042983
#>      rust   4.363978   4.39203   4.621804   4.690928   4.787278   4.916459
#>  neval
#>     10
#>     10
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
#>  ■■■                                6% |  ETA: 24s
#>  ■■■■■■                            18% |  ETA: 21s
#>  ■■■■■■■■■■                        30% |  ETA: 18s
#>  ■■■■■■■■■■■■■■                    42% |  ETA: 15s
#>  ■■■■■■■■■■■■■■■■■                 54% |  ETA: 12s
#>  ■■■■■■■■■■■■■■■■■■■■■             66% |  ETA:  9s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■          78% |  ETA:  6s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% |  ETA:  3s

similarity_matrix <- matrix(
  data = unlist(r_results), 
  nrow = length(random_sets_1), 
  ncol = length(random_sets_2),
  byrow = TRUE
)
tictoc::toc()
#> 25.195 sec elapsed
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
#> 14.108 sec elapsed

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
#> 13.372 sec elapsed

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
#> 8.093 sec elapsed
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
#> 0.048 sec elapsed

all.equal(similarity_matrix, rust_res, tolerance = 1e-15)
#> [1] TRUE
```

Yep, that is MUCH faster… These are some of the functions that are
exposed in `bixverse` and can be integrated into other workflows. There
are many more highly specialised functions that can be used in other
workflows. It’s best to explore the package and the underlying code base
to identify all of these functions.
