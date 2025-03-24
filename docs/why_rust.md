# *Why Rust?*

*Last update: 19.03.2025* </br>

Honest answer? The creation of the package was based on one person wanting to
play around with Rust. To make it useful, the initial playground was stuff that
is relevant for the day-to-day work. Quickly, it became apparent that Rust
allowed for incredible speed gains over R. The results were striking to say the
least.

## *Correlations and co-variance*

```
# Improvements in correlation speeds

set.seed(123L)
no_samples <- 1000
no_features <- 1000
random_data <- matrix(rnorm(no_samples * no_features),
                      nrow = no_samples,
                      ncol = no_features)

# Co-variance
microbenchmark::microbenchmark(
  cov_r <- cov(random_data),
  cov_rust <- rs_covariance(random_data),
  times = 10L
)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r <- cor(random_data),
  cor_res_rust <- rs_cor(random_data, spearman = FALSE),
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
for an M1 Max MacBook Pro (using the Apple Accelerate.framework for R BLAS) for
1000 x 1000 matrices:

**Covariance**

|Implementation|median (ms, 20 runs)|
|--:|--:|
|base R| 515.73|
|Rust|18.66|

*A 27x increase in speed.*

**Pearson correlation**

|Implementation|median (ms, 20 runs)|
|--:|--:|
|base R| 512.83|
|Rust|20.90|

*A 24x increase in speed.*

**Covariance**

|Implementation|median (ms, 20 runs)|
|--:|--:|
|base R|603.00|
|Rust|29.47|

*A 20x increase in speed.*

The speed ups were even more pronounced for larger data matrices...

```
# Larger data matrix
set.seed(123L)
no_samples <- 1000
no_features <- 5000
random_data_2 <- matrix(rnorm(no_samples * no_features),
                        nrow = no_samples,
                        ncol = no_features)

# Pearson correlation
microbenchmark::microbenchmark(
  cor_res_r_2 <- cor(random_data_2),
  cor_res_rust_2 <- rs_cor(random_data_2, spearman = FALSE),
  times = 20L
)
```

**Pearson correlation (on larger matrix)**

|Implementation|median (ms, 20 runs)|
|--:|--:|
|base R|12994.37|
|Rust|419.88|

*In this case, the speed gain is 30x...* </br></br>

Based on these observations, and the reality that matrix algebra is underlying
a lot of computational bioinformatics workflows, the realisation that using
Rust to accelerate computations in R was quite obvious... Also, it's a cool
language and you can be a snob for telling people that you are writing Rust.
Let's be honest...
