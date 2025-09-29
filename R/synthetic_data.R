# synthetic data ---------------------------------------------------------------

## gene expression version -----------------------------------------------------

### simple dge version ---------------------------------------------------------

#' Generates a simple synthetic, pseudo gene expression matrix
#'
#' @description
#' This is a function to create a synthetic, pseudo gene expression matrix with
#' a pre-defined number of total genes, groups and number of genes
#' differentially expressed in the groups. The function will create in a
#' standard setting a 1000 x 90 matrix, and add optionally a small group into
#' the mixture that is slightly more difficult to detect
#'
#' @param no_grps Integer. Number of initial larger groups. Default: 3L.
#' @param per_group Integer. Number of samples per larger groups. Default: 30L.
#' @param total_genes Integer. Number of total genes in the matrix. Default:
#' 1000L.
#' @param no_genes_up Integer. Number of genes per group that are differentially
#' expressed. Default: 100L.
#' @param add_small_group Boolean. Add a smaller group that overlaps with group
#' 1? Default: TRUE.
#' @param size_small_grp Integer. Size of the smaller group.
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `synthetic_matrix_simple` class containing:
#' \itemize{
#'  \item mat - The random matrix
#'  \item diff - List of differentially expressed genes per group
#'  \item group - Look-up vector for sample to group mapping
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' synthetic_GEX <- create_synthetic_signal_matrix()
#' synthetic_signal_mat <- synthetic_GEX$mat
#' }
synthetic_signal_matrix <- function(
  no_grps = 3,
  per_group = 30,
  total_genes = 1000L,
  no_genes_up = 100L,
  add_small_group = TRUE,
  size_small_grp = 5L,
  seed = 10101L
) {
  # hard coded
  mean_initial <- 5
  variance_initial <- 5
  mean_up <- 10
  variance_up <- 5

  # checks
  checkmate::qassert(no_grps, "R1(0,)")
  checkmate::qassert(per_group, "R1(0,)")
  checkmate::qassert(total_genes, "I1(0,)")
  checkmate::qassert(no_genes_up, "I1(0,)")
  checkmate::qassert(add_small_group, "B1")
  checkmate::qassert(size_small_grp, "I1(0,)")
  checkmate::qassert(seed, "I1")

  # matrix generation
  set.seed(seed)
  mat <- matrix(
    rnorm(
      total_genes * per_group * no_grps,
      mean_initial,
      variance_initial
    ),
    nrow = total_genes,
    ncol = per_group * no_grps
  )

  # add row and column names
  rownames(mat) <- paste0("gene", seq_len(total_genes))
  colnames(mat) <- paste0("sample", seq_len(per_group * no_grps))

  # create group factor
  group <- factor(unlist(lapply(seq_len(no_grps), function(x) {
    rep(paste0("group", x), per_group)
  })))
  names(group) <- colnames(mat)

  # create differentially expressed genes
  set.seed(seed + 1)
  diff <- lapply(seq_len(no_grps), function(x) {
    start_idx <- ((x - 1) * no_genes_up) + 1
    end_idx <- ((x - 1) * no_genes_up) + no_genes_up
    diff_genes <- rownames(mat)[start_idx:end_idx]

    # update matrix with differential expression
    group_mask <- group == paste0("group", x)
    mat[diff_genes, group_mask] <<-
      mat[diff_genes, group_mask] +
      rnorm(no_genes_up, mean_up, variance_up)

    return(diff_genes)
  })
  names(diff) <- paste0("group", seq_len(no_grps))

  # add small group if requested
  if (add_small_group) {
    small_group_genes <- (total_genes - no_genes_up):total_genes
    small_group_samples <- seq_len(size_small_grp)

    mat[small_group_genes, small_group_samples] <-
      mat[small_group_genes, small_group_samples] + mean_up

    diff$small_group <- rownames(mat)[small_group_genes]

    group <- as.character(group)
    group[small_group_samples] <- "small_group"
    names(group) <- colnames(mat)
    group <- factor(group)
  }

  # ensure no negative expression values
  mat[mat < 0] <- 0

  # finalise results
  result <- list(
    mat = mat,
    diff = diff,
    group = group
  )

  class(result) <- "synthetic_matrix_simple"

  return(result)
}

### correlation structure version ----------------------------------------------

#' Generates synthetic bulk RNAseq data
#'
#' @description
#' Function generates synthetic bulkRNAseq data with heteroskedasticity (lowly
#' expressed genes show higher variance) and can optionally add gene ~ gene
#' correlations for testing purposes of module detection methods.
#'
#' @param num_samples Integer. Number of samples.
#' @param num_genes Integer. Number of genes.
#' @param add_modules Boolean. Shall gene modules with correlation structures
#' be generated.
#' @param module_sizes Integer vector. Sizes of the different correlation
#' modules. The sum needs be smaller than `num_genes`.
#' @param seed Integer. Seed for reproducibility purposes.
#'
#' @return A `synthetic_bulk_data` class containing:
#' \itemize{
#'  \item counts - The count matrix
#'  \item sparse_counts - A slot for sparse counts that can be added later, see
#'  [bixverse::simulate_dropouts()].
#'  \item module_data - The module membership of the genes
#' }
#'
#' @export
synthetic_bulk_cor_matrix <- function(
  num_samples = 100L,
  num_genes = 1000L,
  add_modules = TRUE,
  module_sizes = c(100L, 100L, 100L),
  seed = 123L
) {
  # checks
  checkmate::qassert(num_samples, "I1")
  checkmate::qassert(num_genes, "I1")
  checkmate::qassert(add_modules, "B1")
  checkmate::qassert(module_sizes, "I+")
  checkmate::qassert(seed, "I1")
  checkmate::assertTRUE(sum(module_sizes) <= num_genes)

  # body
  c(counts, module_membership) %<-%
    rs_generate_bulk_rnaseq(
      num_samples = num_samples,
      num_genes = num_genes,
      seed = seed,
      add_modules = add_modules,
      module_sizes = module_sizes
    )

  rownames(counts) <- sprintf("gene_%i", 1:num_genes)
  colnames(counts) <- sprintf("sample_%i", 1:num_samples)

  module_dt <- data.table::data.table(
    gene = rownames(counts),
    membership = module_membership
  )

  result <- list(counts = counts, sparse_counts = NULL, module_data = module_dt)
  class(result) <- "synthetic_bulk_data"

  return(result)
}

#### sparsification ------------------------------------------------------------

#' Simulate dropouts via different functions on synthetic bulk data
#'
#' @description
#' This function induces expression-level dependent sparsity on the data. The
#' two possible functions are:
#'
#' **Logistic function:**
#'
#' With dropout probability defined as:
#'
#' ```
#' P(dropout) = clamp(1 / (1 + exp(shape * (ln(exp+1) - ln(midpoint+1)))), 0.3, 0.8) * (1 - global_sparsity) + global_sparsity
#' ```
#'
#' with the following characteristics:
#'
#' - Plateaus at global_sparsity dropout for high expression genes
#' - Partial dropout preserves count structure via binomial thinning
#' - Good for preserving variance-mean relationships
#'
#'**Power Decay function:**
#'
#' With dropout probability defined as:
#'
#' ```
#' P(dropout) = (midpoint / (exp + midpoint))^power * scale_factor * (1 - global_sparsity) + global_sparsity
#' ```
#'
#' with the following characteristics:
#'
#' - No plateau - high expression genes get substantial dropout
#' - Complete dropout only (no partial dropout)
#' - More uniform dropout across expression range
#'
#' @param object The `synthetic_bulk_data` class.
#' @param dropout_function String. One of `c("log", "powerdecay")`. These
#' are the two options to cause dropout in the bulk data.
#' @param dropout_midpoint Numeric. Controls the midpoint parameter of the
#' logistic and power decay function.
#' @param dropout_shape Numeric. Controls the shape parameter of the logistic
#' function.
#' @param power_factor Numeric. Controls the power factor of the power decay
#' function.
#' @param global_sparsity Numeric. The global sparsity parameter.
#' @param seed Integer. Random seed for reproducibility purposes.
#'
#' @return `synthetic_bulk_data` with added sparse data.
#'
#' @export
simulate_dropouts <- function(
  object,
  dropout_function = c("log", "powerdecay"),
  dropout_midpoint = 5,
  dropout_shape = 0.5,
  power_factor = 0.3,
  global_sparsity = 0.25,
  seed = 123L
) {
  dropout_function <- match.arg(dropout_function)

  # checks
  checkmate::assertClass(object, "synthetic_bulk_data")
  checkmate::assertChoice(dropout_function, c("log", "powerdecay"))
  checkmate::qassert(dropout_midpoint, "N1")
  checkmate::qassert(dropout_shape, "N1")
  checkmate::qassert(power_factor, "N1[0, 1]")
  checkmate::qassert(seed, "I1")

  sparse_counts <- rs_simulate_dropouts(
    count_mat = object$counts,
    dropout_function = dropout_function,
    dropout_midpoint = dropout_midpoint,
    dropout_shape = dropout_shape,
    power_factor = power_factor,
    global_sparsity = global_sparsity,
    seed = seed
  )

  rownames(sparse_counts) <- rownames(object$counts)
  colnames(sparse_counts) <- colnames(object$counts)

  sparsification_params <- list(
    dropout_function = dropout_function,
    dropout_midpoint = dropout_midpoint,
    dropout_shape = dropout_shape,
    power_factor = power_factor,
    global_sparsity = global_sparsity,
    seed = seed
  )

  object$sparse_counts <- sparse_counts
  attr(object, "dropout_params") <- sparsification_params

  return(object)
}


#' Helper function to calculate the induced sparsity
#'
#' @param object `synthetic_bulk_data` object. You need to have run
#' [bixverse::simulate_dropouts()] for this function to work.
#' @param no_exp_bins Integer. Number of expression bins to check. Defaults to
#' `10L`.
#'
#' @returns A list with various statistics about the sparsity
#' \itemize{
#'  \item original_sparsity - Original proportion of zeroes in the counts.
#'  \item final_sparsity - Sparsity after applying
#'  [bixverse::simulate_dropouts()].
#'  \item added_sparsity - Added sparsity.
#'  \item gene_sparsity_mean - Mean sparsity for the genes.
#'  \item gene_sparsity_sd - SD sparsity form the genes.
#'  \item sample_sparsity_mean - Mean sparsity for the genes.
#'  \item sample_sparsity_sd - SD sparsity for the genes.
#'  \item dropout_by_expression - Dropout per expression bin level.
#' }
#'
#' @export
calculate_sparsity_stats <- function(object, no_exp_bins = 10L) {
  original_zero <- sparse_zero <- NULL

  # checks
  checkmate::assertClass(object, "synthetic_bulk_data")
  checkmate::qassert(no_exp_bins, "I1")

  # early return
  if (is.null(object$sparse_counts)) {
    warning(paste(
      "No sparse counts found.",
      "Did you run simulate_dropouts()? Returning NULL"
    ))
    return(NULL)
  }

  original_mat <- object$counts
  sparse_mat <- object$sparse_counts

  # general sparsity calculations
  c(original_zero, original_row_zero, original_col_zero) %<-%
    rs_count_zeroes(original_mat)

  c(sparse_zero, sparse_row_zero, sparse_col_zero) %<-%
    rs_count_zeroes(sparse_mat)

  n <- length(original_mat)

  original_sparsity <- original_zero / n
  final_sparsity <- sparse_zero / n
  added_sparsity <- final_sparsity - original_sparsity

  # dropout per expression bin
  non_zero_original <- original_mat[original_mat > 0]
  corresponding_dropout <- sparse_mat[original_mat > 0]
  dropout_events <- corresponding_dropout == 0

  expr_bins <- cut(log1p(non_zero_original), breaks = no_exp_bins)
  dropout_by_expr <- tapply(dropout_events, expr_bins, mean)

  stats <- list(
    original_sparsity = original_sparsity,
    final_sparsity = final_sparsity,
    added_sparsity = added_sparsity,
    gene_sparsity_mean = mean(sparse_row_zero),
    gene_sparsity_sd = sd(sparse_row_zero),
    sample_sparsity_mean = mean(sparse_col_zero),
    sample_sparsity_sd = sd(sparse_col_zero),
    dropout_by_expression = dropout_by_expr
  )

  return(stats)
}

## contrastive pca synthetic data ----------------------------------------------

#' @title
#' Generates synthetic data for contrastive PCA exploration.
#'
#' @description
#' Generates three elements: target matrix, background matrix and labels for the
#' target matrix.
#'
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `cpca_synthetic_data` class with the following elements:
#' \itemize{
#'  \item target - The target matrix.
#'  \item background - The background matrix.
#'  \item target_labels - The target labels
#' }
#'
#' @importFrom magrittr `%>%`
#'
#' @export
synthetic_c_pca_data <- function(seed = 10101L) {
  # Checks
  checkmate::qassert(seed, "I1")
  set.seed(seed)

  # Background matrix
  background <- matrix(0, nrow = 400, ncol = 30)
  background[, 1:10] <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)
  background[, 11:20] <- matrix(rnorm(400 * 10, mean = 0, sd = 3), ncol = 10)
  background[, 21:30] <- matrix(rnorm(400 * 10, mean = 0, sd = 1), ncol = 10)

  # Target matrix with 4 clusters
  target <- matrix(0, nrow = 400, ncol = 30)
  target[, 1:10] <- matrix(rnorm(400 * 10, mean = 0, sd = 10), ncol = 10)

  # Group 1
  target[1:100, 11:20] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)
  target[1:100, 21:30] <- matrix(rnorm(100 * 10, mean = 0, sd = 1), ncol = 10)

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

  result <- list(
    target = t(target),
    background = t(background),
    target_labels = target_labels
  )

  class(result) <- "cpca_synthetic_data"

  return(result)
}

## gene module data ------------------------------------------------------------

#' @title
#' Generates synthetic gene module data.
#'
#' @description
#' Generates an artifical matrix of active modules (for a maximum of 8) in
#' samples being active. Also allows for overlap between the modules. Designed
#' to test some of the algorithms.
#'
#' @param n_samples Integer. Number of samples.
#' @param n_genes Integer. Number of genes.
#' @param n_modules Integer. Number of active gene modules. To a maximum of 10.
#' @param overlap_fraction Float. How much the same modules can be active
#' in the same samples.
#' @param seed Integer. Initial random seed for generation of the synthetic
#' data. Default: 10101L.
#'
#' @return A `synthetic_matrix_modules` class with the following items:
#' \itemize{
#'  \item data - The data matrix.
#'  \item metadata - The sample metadata.
#' }
#'
#' @importFrom magrittr `%>%`
#'
#' @export
generate_gene_module_data <- function(
  n_samples = 24L,
  n_genes = 60L,
  n_modules = 4L,
  overlap_fraction = 0.25,
  seed = 10101L
) {
  # checks
  checkmate::qassert(n_samples, "I1")
  checkmate::qassert(n_genes, "I1")
  checkmate::qassert(n_modules, "I1[1, 11)")
  checkmate::qassert(overlap_fraction, "N1(0, 1)")
  checkmate::qassert(seed, "I1")

  # body
  set.seed(seed)

  data <- matrix(0, nrow = n_samples, ncol = n_genes)

  genes_per_module <- n_genes %/% n_modules
  gene_boundaries <- seq(1, n_genes + 1, by = genes_per_module)
  if (length(gene_boundaries) > n_modules + 1) {
    gene_boundaries <- gene_boundaries[1:(n_modules + 1)]
  }
  gene_boundaries[n_modules + 1] <- n_genes + 1 # Ensure last boundary covers all genes

  cells_per_module <- n_samples %/% n_modules
  overlap_size <- round(cells_per_module * overlap_fraction)

  modules <- vector("list", n_modules)
  active_cells <- vector("list", n_modules)

  for (i in 1:n_modules) {
    gene_start <- gene_boundaries[i]
    gene_end <- gene_boundaries[i + 1] - 1

    modules[[i]] <- list(
      start = gene_start,
      end = gene_end
    )

    cell_start <- max(1, (i - 1) * cells_per_module - overlap_size + 1)
    cell_end <- min(n_samples, i * cells_per_module + overlap_size)

    active_cells[[i]] <- cell_start:cell_end
  }

  # Simulate module activities
  for (i in seq_along(modules)) {
    module <- modules[[i]]

    for (cell in active_cells[[i]]) {
      base_activity <- 1.0 + runif(1, -0.2, 0.2)

      for (gene in module$start:module$end) {
        data[cell, gene] <- base_activity + runif(1, -0.1, 0.1)
      }
    }
  }

  # Add background noise
  noise <- matrix(
    runif(n_samples * n_genes, -0.05, 0.05),
    nrow = n_samples,
    ncol = n_genes
  )
  data <- data + noise

  rownames(data) <- sprintf("sample_%i", 1:n_samples)
  colnames(data) <- sprintf("feature_%i", 1:n_genes)

  activity <- purrr::map(active_cells, \(sample_idx) {
    1:n_samples %in% sample_idx
  })

  activity_dt <- data.table::data.table(
    sample_id = sprintf("sample_%i", 1:n_samples)
  )

  for (i in seq_along(activity)) {
    activity_dt[, c(sprintf("module_%i_active", i)) := activity[[i]]]
  }

  result <- list(data = data, meta_data = activity_dt)

  class(result) <- "synthetic_matrix_modules"

  return(result)
}

## single cell -----------------------------------------------------------------

#' Single cell test data
#'
#' @description
#' This function generates synthetic data for single cell test purposes. It has
#' hard-coded parameters and will generate a count matrix of 1000 cells x 100
#' genes, an obs table and a var table. There are three distinct cell types
#' that can be found in the data that each express between 2 to 8 marker genes
#' with higher expression compared to the background. The background genes
#' have some with higher expression and then less and less cells expressing
#' them.
#'
#'
#' @param seed Integer. The seed for the generation of the seed data.
#'
#' @returns List with the following items
#' \itemize{
#'   \item counts - dgRMatrix with cells x genes.
#'   \item obs - data.table that contains the cell information.
#'   \item var - data.table that contains the var information.
#' }
#'
#' @export
generate_single_cell_test_data <- function(seed = 42L) {
  # checks
  checkmate::qassert(seed, "I1")

  # hard-coded parameters
  n_cells = 1000L
  n_genes = 100L
  n_background_genes_exp = c(3L, 15L)
  background_exp_range = c(5L, 10L)

  marker_genes <- list(
    cell_type_1 = list(
      marker_genes = 0:9L,
      marker_exp_range = c(10L, 50L),
      markers_per_cell = c(2L, 8L)
    ),
    cell_type_2 = list(
      marker_genes = 10:19L,
      marker_exp_range = c(10L, 50L),
      markers_per_cell = c(2L, 8L)
    ),
    cell_type_3 = list(
      marker_genes = 20:29L,
      marker_exp_range = c(10L, 50L),
      markers_per_cell = c(2L, 8L)
    )
  )

  data <- rs_synthetic_sc_data_with_cell_types(
    n_cells = n_cells,
    n_genes = n_genes,
    n_background_genes_exp = n_background_genes_exp,
    background_exp_range = n_background_genes_exp,
    cell_configs = marker_genes,
    seed = 42L
  )

  counts <- new(
    "dgRMatrix",
    p = as.integer(data$indptr),
    x = as.numeric(data$data),
    j = as.integer(data$indices),
    Dim = as.integer(c(n_cells, n_genes))
  )

  rownames(counts) <- sprintf("cell_%04d", 1:1000)
  colnames(counts) <- sprintf("gene_%03d", 1:100)

  obs <- data.table(
    cell_id = sprintf("cell_%04d", 1:1000),
    cell_grp = sprintf("cell_type_%i", data$cell_type_indices + 1)
  )

  var <- data.table(
    gene_id = sprintf("gene_%03d", 1:100),
    ensembl_id = sprintf("ens_%03d", 1:100)
  )

  res <- list(
    counts = counts,
    obs = obs,
    var = var
  )
}

# plots ------------------------------------------------------------------------

## simple data -----------------------------------------------------------------

#' Plot the (simple) synthetic gene expression
#'
#' @param x `synthetic_matrix_simple` class. Output from
#' [bixverse::synthetic_signal_matrix()].
#' @param ... Additional params
#'
#' @return A plotted heatmap showing the DEG.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' synthetic_gex <- synthetic_signal_matrix()
#'
#' plot(synthetic_gex)
#' }
plot.synthetic_matrix_simple <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "synthetic_matrix_simple")

  # pull out relevant parts
  mat <- x$mat
  group <- x$group

  # plot
  heatmap(
    mat,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    ColSideColors = rainbow(length(levels(group)))[group],
    labCol = FALSE,
    labRow = FALSE,
    xlab = "Samples",
    ylab = "Features",
    main = "Synthetic gene expression data"
  )
}

## contrastive pca -------------------------------------------------------------

#' @title
#' Plot the contrastive PCA example data
#'
#' @param x `cpca_synthetic_data` class. Output from
#' [bixverse::synthetic_c_pca_data()].
#' @param ... Additional params
#'
#' @return A ggplot showing the two heatmaps from the target and background
#' matrix.
#'
#' @export
#'
#' @import patchwork
#' @import ggplot2
#' @importFrom magrittr `%>%`
plot.cpca_synthetic_data <- function(x, ...) {
  # scope
  . <- NULL

  # checks
  checkmate::assertClass(x, "cpca_synthetic_data")

  # get the data
  target_df <- data.table::as.data.table(
    t(x$target),
    keep.rownames = "samples"
  ) %>%
    .[, grp := x$target_labels] %>%
    data.table::melt(id.vars = c("samples", "grp"))

  background_df <- data.table::as.data.table(
    t(x$background),
    keep.rownames = "samples"
  ) %>%
    data.table::melt(id.vars = "samples")

  # plots
  p1 <- ggplot(
    data = target_df,
    mapping = aes(y = samples, x = variable, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("Target data") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab("Features") +
    ylab("Samples") +
    geom_rug(
      data = unique(target_df[, c("samples", "grp")]),
      mapping = aes(y = samples, color = grp),
      inherit.aes = FALSE,
      sides = "l",
      linewidth = 2,
      length = unit(0.0125, "npc")
    ) +
    labs(fill = "Value:", colour = "Group:")

  p2 <- ggplot(
    data = background_df,
    mapping = aes(y = samples, x = variable, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("Background data") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab("Features") +
    ylab("Samples") +
    labs(fill = "Value:")

  p_final <- p1 +
    p2 +
    patchwork::plot_annotation(
      title = "Contrastive PCA",
      subtitle = "Synthetic data"
    )

  p_final
}
