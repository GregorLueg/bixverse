kline_clustering <- function(data, k = 2, max_iter = 100, tol = 1e-6) {
  # Check if data has exactly 2 columns
  if (ncol(data) != 2) {
    stop("Data must have exactly 2 columns (x and y coordinates)")
  }

  # Convert data to matrix if it's a data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]

  # Initialize cluster assignments randomly
  clusters <- sample(1:k, n, replace = TRUE)

  # Function to compute orthogonal distance from a point to a line
  # Line equation: ax + by + c = 0 (a^2 + b^2 = 1 for normalization)
  # For y = mx + b line: ax + by + c = 0 where a = m, b = -1, c = b
  # We normalize to get a^2 + b^2 = 1
  point_to_line_distance <- function(x, y, m, b) {
    # Convert slope-intercept form to general form (ax + by + c = 0)
    a <- m
    b <- -1
    c <- b

    # Normalize coefficients
    norm_factor <- sqrt(a^2 + b^2)
    a <- a / norm_factor
    b <- b / norm_factor
    c <- c / norm_factor

    # Compute orthogonal distance
    abs(a * x + b * y + c)
  }

  # Initialize variables for iteration
  prev_sse <- Inf
  converged <- FALSE
  iter <- 0

  # Line parameters: slopes and intercepts
  line_params <- matrix(0, nrow = k, ncol = 2)
  colnames(line_params) <- c("slope", "intercept")

  # Iterate until convergence or max iterations
  while (!converged && iter < max_iter) {
    iter <- iter + 1

    # Fit lines to each cluster
    for (i in 1:k) {
      cluster_points <- data[clusters == i, , drop = FALSE]

      # Need at least 2 points to fit a line
      if (nrow(cluster_points) < 2) {
        # If no or few points, initialize with random line
        line_params[i, ] <- c(runif(1, -5, 5), runif(1, -5, 5))
      } else {
        # Fit line using linear regression
        fit <- lm(cluster_points[, 2] ~ cluster_points[, 1])
        line_params[i, ] <- c(coef(fit)[2], coef(fit)[1])
      }
    }

    # Compute distances from each point to each line
    distances <- matrix(0, nrow = n, ncol = k)
    for (i in 1:k) {
      m <- line_params[i, "slope"]
      b <- line_params[i, "intercept"]
      distances[, i] <- point_to_line_distance(x, y, m, b)
    }

    # Assign points to closest line
    new_clusters <- apply(distances, 1, which.min)

    # Compute sum of squared distances (SSE)
    current_sse <- sum(apply(distances, 1, min)^2)

    # Check for convergence
    if (abs(prev_sse - current_sse) < tol) {
      converged <- TRUE
    }

    # Update for next iteration
    prev_sse <- current_sse
    clusters <- new_clusters
  }

  # Compute final distances
  distances <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    m <- line_params[i, "slope"]
    b <- line_params[i, "intercept"]
    distances[, i] <- point_to_line_distance(x, y, m, b)
  }

  final_distances <- numeric(n)
  for (i in 1:n) {
    final_distances[i] <- distances[i, clusters[i]]
  }

  return(list(
    lines = line_params,
    clusters = clusters,
    distances = final_distances,
    iterations = iter,
    converged = converged
  ))
}

# Example usage with simulated data
set.seed(123)

# Generate points along two lines with some noise
generate_line_data <- function(n_points = 100, k = 2, noise = 0.5) {
  # Define k lines
  slopes <- runif(k, -2, 2)
  intercepts <- runif(k, -5, 5)

  x <- numeric(n_points)
  y <- numeric(n_points)
  true_clusters <- numeric(n_points)

  points_per_line <- n_points / k

  for (i in 1:k) {
    idx_start <- round((i - 1) * points_per_line) + 1
    idx_end <- round(i * points_per_line)
    n_line_points <- idx_end - idx_start + 1

    x_range <- runif(n_line_points, -10, 10)
    x[idx_start:idx_end] <- x_range
    y[idx_start:idx_end] <- slopes[i] *
      x_range +
      intercepts[i] +
      rnorm(n_line_points, 0, noise)
    true_clusters[idx_start:idx_end] <- i
  }

  # Shuffle the data
  idx <- sample(n_points)
  x <- x[idx]
  y <- y[idx]
  true_clusters <- true_clusters[idx]

  return(list(
    data = data.frame(x = x, y = y),
    true_clusters = true_clusters,
    true_lines = data.frame(slope = slopes, intercept = intercepts)
  ))
}

# Generate example data
example_data <- generate_line_data(n_points = 200, k = 2, noise = 0.8)

# Run k-line clustering
result <- kline_clustering(example_data$data, k = 2)


as.matrix(ica_stability_res[, c("component_rank", "stability")])


result_ica <- kline_clustering(
  data = as.matrix(ica_stability_res[, c("component_rank", "stability")]),
  k = 2
)


# Plot the results
library(ggplot2)

# Create a data frame for plotting
plot_data <- data.frame(
  x = example_data$data$x,
  y = example_data$data$y,
  cluster = as.factor(result$clusters)
)

ggplot(data = plot_data, mapping = aes(x = x, y = y)) +
  geom_point(aes(col = cluster))

# Generate line points for plotting
line_points <- list()
for (i in 1:2) {
  m <- result$lines[i, "slope"]
  b <- result$lines[i, "intercept"]
  x_range <- seq(min(plot_data$x), max(plot_data$x), length.out = 100)
  y_range <- m * x_range + b
  line_points[[i]] <- data.frame(
    x = x_range,
    y = y_range,
    cluster = as.factor(i)
  )
}
line_df <- do.call(rbind, line_points)

# Create the plot
p <- ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
  geom_point(alpha = 0.7) +
  geom_line(data = line_df, aes(x = x, y = y, color = cluster), size = 1) +
  theme_minimal() +
  ggtitle("K-Line Clustering Results (k = 2)") +
  labs(x = "X", y = "Y", color = "Cluster")

print(p)
