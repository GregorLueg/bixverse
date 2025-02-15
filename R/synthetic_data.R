# co-expression synthetic data ----

## cPCA synthetic data ----

#' Generates synthetic data for cPCA exploration.
#'
#' @description Generates three elements: target matrix, background matrix
#' and labels for the target matrix.
#'
#' @param seed Integer. Initial random seed for generation of the synthetic data. Default: 10101L.
#'
#' @return A list. First elements is the target matrix, second element the background matrix and the
#' third one being a factor with the cluster labels of the target matrix.
#' @importFrom magrittr `%>%`
#'
#' @export
create_synthetic_cPCA_data <- function(seed = 10101L) {
  # Checks
  checkmate::qassert(seed, "I1")
  set.seed(seed)

  # Background matrix
  background <- matrix(0, nrow = 400, ncol = 30)
  background[, 1:10]  <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  background[, 11:20] <- matrix(rnorm(400 * 10, mean = 0, sd = 3), ncol = 10)
  background[, 21:30] <- matrix(rnorm(400 * 10, mean = 0, sd = 1), ncol = 10)

  # Target matrix with 4 clusters
  target <- matrix(0, nrow = 400, ncol = 30)
  target[, 1:10]  <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)

  # Group 1
  target[1:100, 11:20]  <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[1:100, 21:30]  <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

  # Group 2
  target[101:200, 11:20] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[101:200, 21:30] <- matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)

  # Group 3
  target[201:300, 11:20] <- matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[201:300, 21:30] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

  # Group 4
  target[301:400, 11:20] <- matrix(rnorm(100 * 10, mean = 6, sd = 1), ncol = 10)
  target[301:400, 21:30] <- matrix(rnorm(100 * 10, mean = 3, sd = 1), ncol = 10)

  # Labels
  target_labels <- factor(as.character(rep(0:3, each = 100)))

  # Assign row and column names
  rownames(background) <- rownames(target) <- paste0("sample_", 1:400)
  colnames(background) <- colnames(target) <- paste0("feature_", 1:30)

  list(
    target = t(target),
    background = t(background),
    target_labels = target_labels
  )
}
