library(magrittr)

## Synthetic data ----

create_synthetic_cPCA_data <- function(seed = 10101L) {
  # Checks
  checkmate::qassert(seed, "I1")
  set.seed(seed)
  # Background matrix
  background <- matrix(0, nrow = 400, ncol = 30)
  background[, 1:10] <-
    matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  background[, 11:20] <-
    matrix(rnorm(400 * 10, mean = 0, sd = 3), ncol = 10)
  background[, 21:30] <-
    matrix(rnorm(400 * 10, mean = 0, sd = 1), ncol = 10)
  # Target matrix with 4 clusters
  # The first 10 features have the same variance as the background data set
  # And generate strong background variance, hiding the actual signal
  target <- matrix(0, nrow = 400, ncol = 30)
  target[, 1:10] <-
    matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  # Group 1
  target[1:100, 11:20] <-
    matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[1:100, 21:30] <-
    matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  # Group 2
  target[101:200, 11:20] <-
    matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[101:200, 21:30] <-
    matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)
  # Group 3
  target[201:300, 11:20] <-
    matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[201:300, 21:30] <-
    matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  # Group 4
  target[301:400, 11:20] <-
    matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[301:400, 21:30] <-
    matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)
  # Labels
  target_labels <-
    c(rep(0, 100), rep(1, 100), rep(2, 100), rep(3, 100)) %>%
    as.character() %>%
    as.factor()
  rownames(background) <- rownames(target) <- paste("sample", 1:400, sep = "_")
  colnames(background) <- colnames(target) <- paste("feature", 1:30, sep = "_")
  # Return the results
  res <- list(
    target = t(target),
    background = t(background),
    target_labels = target_labels
  )
}

## Generate the data ----

cPCA_data = create_synthetic_cPCA_data()




