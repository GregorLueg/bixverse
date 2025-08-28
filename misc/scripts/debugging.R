# functions --------------------------------------------------------------------

generate_bulk_rnaseq <- function(
  num_samples = 20L,
  num_genes = 10000L,
  seed = 123L,
  add_modules = FALSE,
  module_sizes = c(300, 250, 200, 300, 500)
) {
  # Checks
  checkmate::assertInteger(num_samples)
  checkmate::assertInteger(num_genes)
  checkmate::assertInteger(seed)
  checkmate::assertLogical(add_modules)
  checkmate::assertVector(module_sizes)
  checkmate::assert(sum(module_sizes) <= num_genes)

  # Set seed for reproducibility
  set.seed(seed)

  # Generate mean expression values for each gene (typically follows a gamma distribution)
  mean_expression <- rgamma(num_genes, shape = 5, scale = 10)

  # Generate dispersion parameters (typically inversely related to mean)
  # Higher dispersion = more variance
  dispersion <- 1 / (0.5 + mean_expression)

  # Initialize matrix to store counts
  count_matrix <- matrix(0, nrow = num_genes, ncol = num_samples)
  rownames(count_matrix) <- paste0("gene", 1:num_genes)
  colnames(count_matrix) <- paste0("sample", 1:num_samples)

  # Fill the matrix with negative binomial distributed counts
  for (i in 1:num_genes) {
    count_matrix[i, ] <- rnbinom(
      num_samples,
      size = 1 / dispersion[i],
      mu = mean_expression[i]
    )
  }

  # Create module assignments for each gene
  module_assignment <- rep(0, num_genes) # 0 represents background genes

  if (add_modules) {
    # Remaining background genes
    num_background_genes <- num_genes - sum(module_sizes)

    # Create module assignments for each gene
    gene_idx <- 1
    for (module in 1:length(module_sizes)) {
      module_end <- gene_idx + module_sizes[module] - 1
      module_assignment[gene_idx:module_end] <- module
      gene_idx <- module_end + 1
    }

    # Create latent factors for each module (these represent shared biology)
    module_factors <- matrix(
      0,
      nrow = length(module_sizes),
      ncol = num_samples
    )
    for (i in 1:length(module_sizes)) {
      # Generate random factor with some variation between samples
      module_factors[i, ] <- rnorm(num_samples, mean = 1, sd = 0.7)
    }

    # Generate expression data that includes both module-specific correlation and individual variation
    for (i in 1:num_genes) {
      # For genes in modules, add correlated component
      if (module_assignment[i] > 0) {
        module_id <- module_assignment[i]

        # Generate correlation strength (different for each gene in module)
        # Higher values = stronger correlation with module
        correlation_strength <- rbeta(1, 5, 2) # Most genes have strong correlation

        # Create signal by combining independent component with module factor
        module_signal <- exp(
          correlation_strength * module_factors[module_id, ]
        )

        # Scale signal to maintain approximately the same mean
        module_signal <- module_signal *
          (mean_expression[i] / mean(module_signal))

        # Generate counts incorporating the module effect
        count_matrix[i, ] <- rnbinom(
          num_samples,
          size = 1 / dispersion[i],
          mu = module_signal
        )
      } else {
        # Background genes just get independent negative binomial counts
        count_matrix[i, ] <- count_matrix[i, ]
      }
    }
    names(module_assignment) <- paste0("gene", 1:num_genes)
  }
  gene_info <- data.table(
    gene = rownames(count_matrix),
    module = module_assignment
  )
  attr(count_matrix, "gene_info") <- gene_info
  attr(count_matrix, "num_genes") <- num_genes
  attr(count_matrix, "num_samples") <- num_samples
  return(count_matrix)
}

simulate_dropouts <- function(
  count_matrix,
  dropout_midpoint = 1,
  dropout_shape = 1,
  global_sparsity = 0
) {
  # Ensure matrix is numeric
  bulk_matrix <- as.matrix(count_matrix)

  # Get dimensions
  n_genes <- nrow(bulk_matrix)
  n_cells <- ncol(bulk_matrix)

  # Initialize output matrix
  dropout_matrix <- bulk_matrix
  gene_info <- attributes(bulk_matrix)$gene_info

  # Apply dropouts
  for (i in 1:n_genes) {
    for (j in 1:n_cells) {
      expr_value <- bulk_matrix[i, j]

      # Calculate dropout probability for this expression level
      dropout_prob <- calculate_dropout_prob(
        expr_value,
        dropout_midpoint,
        dropout_shape,
        global_sparsity
      )

      # Randomly determine if this value drops out
      if (runif(1) < dropout_prob) {
        dropout_matrix[i, j] <- 0
      }
    }
  }

  attr(dropout_matrix, "gene_info") <- attributes(count_matrix)$gene_info
  attr(dropout_matrix, "parameters") <- list(
    global_sparsity = global_sparsity,
    dropout_midpoint = dropout_midpoint,
    dropout_shape = dropout_shape
  )
  return(dropout_matrix)
}

calculate_dropout_prob <- function(
  expression,
  midpoint,
  shape,
  global_sparsity
) {
  # Logistic function: P(dropout) = 1 / (1 + exp(shape * (expression - midpoint)))
  # This gives high dropout probability for low expression
  prob <- 1 / (1 + exp(shape * (log1p(expression) - log1p(midpoint))))

  # Apply global sparsity factor
  # This increases dropout probability across all expression levels
  prob <- prob + global_sparsity * (1 - prob)

  return(pmin(prob, 1)) # Ensure probability doesn't exceed 1
}

# test bulk rnaseq data --------------------------------------------------------

tictoc::tic()
test = rs_generate_bulk_rnaseq(
  num_samples = 500L,
  num_genes = 10000L,
  seed = 123L,
  add_modules = TRUE,
  module_sizes = NULL
)
tictoc::toc()

tictoc::tic()
dropout_data <- rs_simulate_dropouts(
  count_mat = test$counts,
  dropout_midpoint = 10,
  dropout_shape = 1.5,
  global_sparsity = 0.1,
  seed = 123L
)
tictoc::toc()

norm_counts = edgeR::cpm(test$counts, log = TRUE)

tictoc::tic()
test2 = generate_bulk_rnaseq(
  num_samples = 500L,
  num_genes = 10000L,
  add_modules = TRUE
)

dropout_data_2 <- simulate_dropouts(
  test2,
  dropout_midpoint = 10,
  dropout_shape = 1.5,
  global_sparsity = 0.1
)
tictoc::toc()

plot(rowMeans(norm_counts_2), matrixStats::rowSds(norm_counts_2))
