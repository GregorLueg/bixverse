# *Why Rust?*

*Last update: 19.03.2025* </br>

Honest answer? The creation of the package was based on one person wanting to
play around with Rust. To make it useful, the initial playground was stuff that
is relevant for the day-to-day work. Quickly, it became apparent that Rust
allowed for incredible speed gains over R. The results were striking...

## *Correlations and co-variance*

```
# Improvements in correlation speeds

set.seed(123L)
no_samples <- 1000
no_features <- 1000
random_data <- matrix(rnorm(no_samples * no_features),
                      nrow = no_samples,
                      ncol = no_features)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data),
  cor_res_rust <- rs_cor(random_data, spearman = FALSE),
  times = 10L
)

# Co-variance
microbenchmark::microbenchmark(
  cov_r <- cov(random_data),
  cov_rust <- rs_covariance(random_data),
  times = 10L
)

# Spearman correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data, method = 'spearman'),
  cor_res_rust <- rs_cor(random_data, spearman = TRUE),
  times = 10L
)
```

In this simple example, the *rs_* variants had the following increases in speed
for an M1 Max MacBook Pro (using the Apple Accelerate.framework for R BLAS):

**Pearson correlation**

|Implementation|mean (ms, 10 runs)|
|:--|--:|
|base R| 432.41|
|Rust|27.39|

*A 15x increase in speed.*

**Covariance**

|Implementation|mean (ms, 10 runs)|
|:--|--:|
|base R| 435.66|
|Rust|18.26|

*A 20x increase in speed.*

**Covariance**

|Implementation|mean (ms, 10 runs)|
|:--|--:|
|base R| 534.10|
|Rust|25.38|

*A 20x increase in speed.*

## *Other examples* 

To be added...

